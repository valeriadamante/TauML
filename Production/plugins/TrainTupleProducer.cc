/*! Creates tuple for tau analysis.
*/

#include "Compression.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"

#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "TauML/Analysis/include/TrainTuple.h"
#include "TauML/Analysis/include/SummaryTuple.h"
#include "TauML/Analysis/include/TauIdResults.h"
#include "TauML/Production/include/GenTruthTools.h"
#include "TauML/Production/include/TauAnalysis.h"
#include "TauML/Production/include/MuonHitMatch.h"
#include "TauML/Production/include/TauJet.h"

#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterAssociation.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameterFwd.h"
#include "DataFormats/TauReco/interface/PFTauTransverseImpactParameter.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
// #include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

namespace tau_analysis {

struct TrainTupleProducerData {
    using clock = std::chrono::system_clock;

    const clock::time_point start;
    train_tuple::TrainTuple trainTuple;
    tau_tuple::SummaryTuple summaryTuple;
    std::mutex mutex;

private:
    size_t n_producers;

    TrainTupleProducerData(TFile& file) :
        start(clock::now()),
        trainTuple("taus", &file, false),
        summaryTuple("summary", &file, false),
        n_producers(0)
    {
        summaryTuple().numberOfProcessedEvents = 0;
    }

    ~TrainTupleProducerData() {}

public:

    static TrainTupleProducerData* RequestGlobalData()
    {
        TrainTupleProducerData* data = GetGlobalData();
        if(data == nullptr)
            throw cms::Exception("TrainTupleProducerData") << "Request after all data copies were released.";
        {
            std::lock_guard<std::mutex> lock(data->mutex);
            ++data->n_producers;
            std::cout << "New request of TrainTupleProducerData. Total number of producers = " << data->n_producers
                      << "." << std::endl;
        }
        return data;
    }

    static void ReleaseGlobalData()
    {
        TrainTupleProducerData*& data = GetGlobalData();
        if(data == nullptr)
            throw cms::Exception("TrainTupleProducerData") << "Another release after all data copies were released.";
        {
            std::lock_guard<std::mutex> lock(data->mutex);
            if(!data->n_producers)
                throw cms::Exception("TrainTupleProducerData") << "Release before any request.";
            --data->n_producers;
            std::cout << "TrainTupleProducerData has been released. Total number of producers = " << data->n_producers
                      << "." << std::endl;
            if(!data->n_producers) {
                data->trainTuple.Write();
                const auto stop = clock::now();
                data->summaryTuple().exeTime = static_cast<unsigned>(
                            std::chrono::duration_cast<std::chrono::seconds>(stop - data->start).count());
                data->summaryTuple.Fill();
                data->summaryTuple.Write();
                delete data;
                data = nullptr;
                std::cout << "TrainTupleProducerData has been destroyed." << std::endl;
            }
        }

    }

private:
    static TrainTupleProducerData*& GetGlobalData()
    {
        static TrainTupleProducerData* data = InitializeGlobalData();
        return data;
    }

    static TrainTupleProducerData* InitializeGlobalData()
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kZLIB);
        file.SetCompressionLevel(9);
        TrainTupleProducerData* data = new TrainTupleProducerData(file);
        std::cout << "TrainTupleProducerData has been created." << std::endl;
        return data;
    }
};

class TrainTupleProducer : public edm::EDAnalyzer {
public:
    using TauDiscriminator = reco::PFTauDiscriminator;

    TrainTupleProducer(const edm::ParameterSet& cfg) :
        isMC(cfg.getParameter<bool>("isMC")),
        storeJetsWithoutTau(cfg.getParameter<bool>("storeJetsWithoutTau")),
        requireGenMatch(cfg.getParameter<bool>("requireGenMatch")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(consumes<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        rho_token(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
        electrons_token(consumes<std::vector<reco::RecoEcalCandidate>>(cfg.getParameter<edm::InputTag>("electrons"))),
        muons_token(consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        taus_token(consumes<std::vector<reco::PFTau>>(cfg.getParameter<edm::InputTag>("taus"))),
        l1Taus_token(consumes<l1t::TauBxCollection>(cfg.getParameter<edm::InputTag>("l1taus"))),
        caloTowers_token(consumes<CaloTowerCollection>(cfg.getParameter<edm::InputTag>("caloTowers"))),
        pixelTracks_token(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("pixelTracks"))),
        caloTaus_token(consumes<reco::CaloJetCollection>(cfg.getParameter<edm::InputTag>("caloTaus"))),
        cands_token(consumes<std::vector<reco::PFCandidate>>(cfg.getParameter<edm::InputTag>("pfCandidates"))),
        decayMode_token(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("decayModeFindingNewDM"))),
        chargedIsoPtSum_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("chargedIsoPtSum"))),
        chargedIsoPtSumdR03_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("chargedIsoPtSumdR03"))),
        neutralIsoPtSum_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("neutralIsoPtSum"))),
        neutralIsoPtSumdR03_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("neutralIsoPtSumdR03"))),
        footprintCorrection_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("footprintCorrection"))),
        footprintCorrectiondR03_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("footprintCorrectiondR03"))),
        neutralIsoPtSumWeight_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("neutralIsoPtSumWeight"))),
        neutralIsoPtSumWeightdR03_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("neutralIsoPtSumWeightdR03"))),
        photonPtSumOutsideSignalCone_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("photonPtSumOutsideSignalCone"))),
        photonPtSumOutsideSignalConedR03_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("photonPtSumOutsideSignalConedR03"))),
        puCorrPtSum_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("puCorrPtSum"))),
        PFTauTransverseImpactParameters_token(consumes<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>>(cfg.getParameter<edm::InputTag>("pfTauTransverseImpactParameters"))),
        deepTauVSe_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("deepTauVSe"))),
        deepTauVSmu_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("deepTauVSmu"))),
        deepTauVSjet_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("deepTauVSjet"))),
        looseIsoAbs_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("looseIsoAbs"))),
        looseIsoRel_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("looseIsoRel"))),
        mediumIsoAbs_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("mediumIsoAbs"))),
        mediumIsoRel_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("mediumIsoRel"))),
        tightIsoAbs_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("tightIsoAbs"))),
        tightIsoRel_inputToken(consumes<TauDiscriminator>(cfg.getParameter<edm::InputTag>("tightIsoRel"))),
        data(TrainTupleProducerData::RequestGlobalData()),
        trainTuple(data->trainTuple),
        summaryTuple(data->summaryTuple)
    {
        builderSetup.minJetPt = cfg.getParameter<double>("minJetPt");
        builderSetup.maxJetEta = cfg.getParameter<double>("maxJetEta");
        builderSetup.forceTauJetMatch = cfg.getParameter<bool>("forceTauJetMatch");
        builderSetup.useOnlyTauObjectMatch = !storeJetsWithoutTau;
        builderSetup.tauJetMatchDeltaR2Threshold = std::pow(cfg.getParameter<double>("tauJetMatchDeltaRThreshold"), 2);
        builderSetup.objectMatchDeltaR2ThresholdJet =
                std::pow(cfg.getParameter<double>("objectMatchDeltaRThresholdJet"), 2);
        builderSetup.objectMatchDeltaR2ThresholdTau =
                std::pow(cfg.getParameter<double>("objectMatchDeltaRThresholdTau"), 2);
    }

private:
    static constexpr float default_value = train_tuple::DefaultFillValue<float>();
    static constexpr int default_int_value = train_tuple::DefaultFillValue<int>();
    static constexpr int default_unsigned_value = train_tuple::DefaultFillValue<unsigned>();

    virtual void analyze(const edm::Event& event, const edm::EventSetup&) override
    {
        std::lock_guard<std::mutex> lock(data->mutex);
        summaryTuple().numberOfProcessedEvents++;

        trainTuple().run  = event.id().run();
        trainTuple().lumi = event.id().luminosityBlock();
        trainTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::Vertex>> vertices;
        event.getByToken(vertices_token, vertices);
        trainTuple().npv = static_cast<int>(vertices->size());
        edm::Handle<double> rho;
        event.getByToken(rho_token, rho);
        trainTuple().rho = static_cast<float>(*rho);

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            trainTuple().genEventWeight = static_cast<float>(genEvent->weight());

            edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
            event.getByToken(puInfo_token, puInfo);
            trainTuple().npu = analysis::gen_truth::GetNumberOfPileUpInteractions(puInfo);
        }

        const auto& PV = vertices->at(0);
        trainTuple().pv_x = static_cast<float>(PV.position().x());
        trainTuple().pv_y = static_cast<float>(PV.position().y());
        trainTuple().pv_z = static_cast<float>(PV.position().z());
        trainTuple().pv_chi2 = static_cast<float>(PV.chi2());
        trainTuple().pv_ndof = static_cast<float>(PV.ndof());

        edm::Handle<std::vector<reco::RecoEcalCandidate>> electrons;
        event.getByToken(electrons_token, electrons);

        edm::Handle<reco::MuonCollection> muons;
        event.getByToken(muons_token, muons);

        edm::Handle<std::vector<reco::PFTau>> taus;
        event.getByToken(taus_token, taus);

        edm::Handle<l1t::TauBxCollection> l1Taus;
        event.getByToken(l1Taus_token, l1Taus);

        edm::Handle<CaloTowerCollection> caloTowers;
        event.getByToken(caloTowers_token, caloTowers);

        edm::Handle<reco::TrackCollection> pixelTracks;
        event.getByToken(pixelTracks_token, pixelTracks);

        edm::Handle<reco::CaloJetCollection> caloTaus;
        event.getByToken(caloTaus_token, caloTaus);
        pat::JetCollection jets;

        edm::Handle<std::vector<reco::PFCandidate>> cands;
        event.getByToken(cands_token, cands);

        edm::Handle<reco::PFTauDiscriminator> decayModesNew;
        event.getByToken(decayMode_token, decayModesNew);

        edm::Handle<TauDiscriminator> chargedIsoPtSum;
        event.getByToken(chargedIsoPtSum_inputToken, chargedIsoPtSum);

        edm::Handle<TauDiscriminator> chargedIsoPtSumdR03;
        event.getByToken(chargedIsoPtSumdR03_inputToken, chargedIsoPtSumdR03);

        edm::Handle<TauDiscriminator> neutralIsoPtSum;
        event.getByToken(neutralIsoPtSum_inputToken, neutralIsoPtSum);

        edm::Handle<TauDiscriminator> neutralIsoPtSumdR03;
        event.getByToken(neutralIsoPtSumdR03_inputToken, neutralIsoPtSumdR03);

        edm::Handle<TauDiscriminator> footprintCorrection;
        event.getByToken(footprintCorrection_inputToken, footprintCorrection);

        edm::Handle<TauDiscriminator> footprintCorrectiondR03;
        event.getByToken(footprintCorrectiondR03_inputToken, footprintCorrectiondR03);

        edm::Handle<TauDiscriminator> neutralIsoPtSumWeight;
        event.getByToken(neutralIsoPtSumWeight_inputToken, neutralIsoPtSumWeight);

        edm::Handle<TauDiscriminator> neutralIsoPtSumWeightdR03;
        event.getByToken(neutralIsoPtSumWeightdR03_inputToken, neutralIsoPtSumWeightdR03);

        edm::Handle<TauDiscriminator> photonPtSumOutsideSignalCone;
        event.getByToken(photonPtSumOutsideSignalCone_inputToken, photonPtSumOutsideSignalCone);

        edm::Handle<TauDiscriminator> photonPtSumOutsideSignalConedR03;
        event.getByToken(photonPtSumOutsideSignalConedR03_inputToken, photonPtSumOutsideSignalConedR03);

        edm::Handle<TauDiscriminator> puCorrPtSum;
        event.getByToken(puCorrPtSum_inputToken, puCorrPtSum);

        edm::Handle<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>> PFTauTransverseImpactParameters;
        event.getByToken(PFTauTransverseImpactParameters_token, PFTauTransverseImpactParameters);

        edm::Handle<TauDiscriminator> looseIsoAbs;
        event.getByToken(looseIsoAbs_inputToken, looseIsoAbs);

        edm::Handle<TauDiscriminator> looseIsoRel;
        event.getByToken(looseIsoRel_inputToken, looseIsoRel);

        edm::Handle<TauDiscriminator> mediumIsoAbs;
        event.getByToken(mediumIsoAbs_inputToken, mediumIsoAbs);

        edm::Handle<TauDiscriminator> mediumIsoRel;
        event.getByToken(mediumIsoRel_inputToken, mediumIsoRel);

        edm::Handle<TauDiscriminator> tightIsoAbs;
        event.getByToken(tightIsoAbs_inputToken, tightIsoAbs);

        edm::Handle<TauDiscriminator> tightIsoRel;
        event.getByToken(tightIsoRel_inputToken, tightIsoRel);

        edm::Handle<TauDiscriminator> deepTau_VSe;
        event.getByToken(deepTauVSe_inputToken, deepTau_VSe);

        edm::Handle<TauDiscriminator> deepTau_VSmu;
        event.getByToken(deepTauVSmu_inputToken, deepTau_VSmu);

        edm::Handle<TauDiscriminator> deepTau_VSjet;
        event.getByToken(deepTauVSjet_inputToken, deepTau_VSjet);

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
        if(isMC)
            event.getByToken(genParticles_token, hGenParticles);

        auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

        //check all inputs!!
        for(unsigned n = 0; n < taus->size(); ++n) {
            const reco::PFTau& tau = taus->at(n);
            if(!(tau.polarP4().pt() > 20 && std::abs(tau.polarP4().eta()) < 2.3)) continue;

            const reco::PFCandidate* leadChargedHadrCand = dynamic_cast<const reco::PFCandidate*>(tau.leadChargedHadrCand().get());

            trainTuple().tau_index.push_back(n);
            trainTuple().tau_pt.push_back(static_cast<float>(tau.polarP4().pt()));
            trainTuple().tau_eta.push_back(static_cast<float>(tau.polarP4().eta()));
            trainTuple().tau_phi.push_back(static_cast<float>(tau.polarP4().phi()));
            trainTuple().tau_mass.push_back(static_cast<float>(tau.polarP4().mass()));
            trainTuple().tau_charge.push_back(tau.charge());



            trainTuple().tau_decayMode.push_back(tau.decayMode());
            trainTuple().tau_decayModeFindingNewDMs.push_back(decayModesNew->value(n));
            trainTuple().tau_chargedIsoPtSum.push_back(chargedIsoPtSum->value(n));
            trainTuple().tau_chargedIsoPtSumdR03.push_back(chargedIsoPtSumdR03->value(n));
            trainTuple().tau_footprintCorrection.push_back(footprintCorrection->value(n));
            trainTuple().tau_footprintCorrectiondR03.push_back(footprintCorrectiondR03->value(n));
            trainTuple().tau_neutralIsoPtSum.push_back(neutralIsoPtSum->value(n));
            trainTuple().tau_neutralIsoPtSumWeight.push_back(neutralIsoPtSumWeight->value(n));
            trainTuple().tau_neutralIsoPtSumWeightdR03.push_back(neutralIsoPtSumWeightdR03->value(n));
            trainTuple().tau_neutralIsoPtSumdR03.push_back(neutralIsoPtSumdR03->value(n));
            trainTuple().tau_photonPtSumOutsideSignalCone.push_back(photonPtSumOutsideSignalCone->value(n));
            trainTuple().tau_photonPtSumOutsideSignalConedR03.push_back(photonPtSumOutsideSignalConedR03->value(n));
            trainTuple().tau_puCorrPtSum.push_back(puCorrPtSum->value(n));

            analysis::TauIdResults cutBasedRelIso;
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Loose,looseIsoRel->value(n) > 0.5);
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Medium,mediumIsoRel->value(n) > 0.5);
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Tight,tightIsoRel->value(n) > 0.5);
            trainTuple().tau_cutBasedRelIso.push_back(cutBasedRelIso.GetResultBits());

            analysis::TauIdResults cutBasedAbsIso;
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Loose,looseIsoAbs->value(n) > 0.5);
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Medium,mediumIsoAbs->value(n) > 0.5);
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Tight,tightIsoAbs->value(n) > 0.5);
            trainTuple().tau_cutBasedAbsIso.push_back(cutBasedAbsIso.GetResultBits());

            trainTuple().tau_byDeepTau2017v2VSeraw.push_back(static_cast<float>(deepTau_VSe->value(n)));
            trainTuple().tau_byDeepTau2017v2VSmuraw.push_back(static_cast<float>(deepTau_VSmu->value(n)));
            trainTuple().tau_byDeepTau2017v2VSjetraw.push_back(static_cast<float>(deepTau_VSjet->value(n)));

            auto impactParam = PFTauTransverseImpactParameters->value(n);

            trainTuple().tau_dxy.push_back(impactParam->dxy());
            trainTuple().tau_dxy_error.push_back(impactParam->dxy_error());
            trainTuple().tau_ip3d.push_back(impactParam->ip3d());
            trainTuple().tau_ip3d_error.push_back(impactParam->ip3d_error());
            const bool has_sv = impactParam->hasSecondaryVertex();
            trainTuple().tau_hasSecondaryVertex.push_back(impactParam->hasSecondaryVertex());
            trainTuple().tau_sv_x.push_back(has_sv ? impactParam->secondaryVertexPos().x() : default_value);
            trainTuple().tau_sv_y.push_back(has_sv ? impactParam->secondaryVertexPos().y() : default_value);
            trainTuple().tau_sv_z.push_back(has_sv ? impactParam->secondaryVertexPos().z() : default_value);
            trainTuple().tau_flightLength_x.push_back(impactParam->flightLength().x());
            trainTuple().tau_flightLength_y.push_back(impactParam->flightLength().y());
            trainTuple().tau_flightLength_z.push_back(impactParam->flightLength().z());
            trainTuple().tau_flightLength_sig.push_back(impactParam->flightLengthSig());


            trainTuple().tau_dz.push_back(leadChargedHadrCand  && leadChargedHadrCand->bestTrack() != nullptr ? leadChargedHadrCand->bestTrack()->dz() : default_value);
            trainTuple().tau_dz_error.push_back(leadChargedHadrCand ? leadChargedHadrCand->dzError() : default_value);

            trainTuple().tau_pt_weighted_deta_strip.push_back(reco::tau::pt_weighted_deta_strip(tau, tau.decayMode()));
            trainTuple().tau_pt_weighted_dphi_strip.push_back(reco::tau::pt_weighted_dphi_strip(tau, tau.decayMode()));
            trainTuple().tau_pt_weighted_dr_signal.push_back(reco::tau::pt_weighted_dr_signal(tau, tau.decayMode()));
            trainTuple().tau_pt_weighted_dr_iso.push_back(reco::tau::pt_weighted_dr_iso(tau, tau.decayMode()));
            trainTuple().tau_leadingTrackNormChi2.push_back(reco::tau::lead_track_chi2(tau));
            trainTuple().tau_e_ratio.push_back(reco::tau::eratio(tau));
            trainTuple().tau_gj_angle_diff.push_back(CalculateGottfriedJacksonAngleDifference(tau));
            trainTuple().tau_n_photons.push_back(static_cast<int>(reco::tau::n_photons_total(tau)));

            float emFraction = -1.;
            float myHCALenergy = 0.;
            float myECALenergy = 0.;
            if(leadChargedHadrCand && leadChargedHadrCand->bestTrack() != nullptr){
                for (const auto& isoPFCand : tau.isolationPFCands()) {
                      myHCALenergy += isoPFCand->hcalEnergy();
                      myECALenergy += isoPFCand->ecalEnergy();
                }
                for (const auto& signalPFCand : tau.signalPFCands()) {
                    myHCALenergy += signalPFCand->hcalEnergy();
                    myECALenergy += signalPFCand->ecalEnergy();
                }
                if (myHCALenergy + myECALenergy != 0.) {
                  emFraction = myECALenergy / (myHCALenergy + myECALenergy);
                }
            }
            trainTuple().tau_emFraction.push_back(emFraction);
            trainTuple().tau_inside_ecal_crack.push_back(IsInEcalCrack(tau.p4().Eta()));

            const std::vector<reco::CandidatePtr>& signalCands = tau.signalCands();
            float leadChargedCandPt = -99;
            float leadChargedCandEtaAtEcalEntrance = -99;
            for (const auto& it : signalCands) {
                const reco::PFCandidate* icand = dynamic_cast<const reco::PFCandidate*>(it.get());
                if (icand != nullptr) {
                      const reco::Track* track = nullptr;
                      if (icand->trackRef().isNonnull())
                        track = icand->trackRef().get();
                      else if (icand->muonRef().isNonnull() && icand->muonRef()->innerTrack().isNonnull())
                        track = icand->muonRef()->innerTrack().get();
                      else if (icand->muonRef().isNonnull() && icand->muonRef()->globalTrack().isNonnull())
                        track = icand->muonRef()->globalTrack().get();
                      else if (icand->muonRef().isNonnull() && icand->muonRef()->outerTrack().isNonnull())
                        track = icand->muonRef()->outerTrack().get();
                      else if (icand->gsfTrackRef().isNonnull())
                        track = icand->gsfTrackRef().get();
                      if (track) {
                        if (track->pt() > leadChargedCandPt) {
                          leadChargedCandEtaAtEcalEntrance = icand->positionAtECALEntrance().eta();
                          leadChargedCandPt = track->pt();
                        }
                      }
                }
            }
            trainTuple().tau_leadChargedCand_etaAtEcalEntrance.push_back(leadChargedCandEtaAtEcalEntrance);


        }

        //functions
        FillGenMatchResult(*genParticles);

        //check inputs
        FillL1Objects(*l1Taus);
        FillPixelTracks(*pixelTracks);
        FillCaloTowers(*caloTowers);
        FillCaloTaus(*caloTaus);
        //change input
        FillPFCandidates(*cands);
        FillElectrons(*electrons);
        FillMuons(*muons,*cands);

        trainTuple.Fill();
    }

    virtual void endJob() override
    {
        TrainTupleProducerData::ReleaseGlobalData();
    }

private:



    void FillGenMatchResult(const std::vector<reco::GenParticle>& genParticles)
    {
        std::vector<analysis::gen_truth::LeptonMatchResult> lepton_results = analysis::gen_truth::CollectGenLeptons(genParticles);
        for(unsigned n = 0; n < lepton_results.size(); ++n){
            analysis::gen_truth::LeptonMatchResult leptonMatch = lepton_results.at(n);
            trainTuple().lepton_gen_match.push_back(static_cast<int>(leptonMatch.match));
            trainTuple().lepton_genLast_charge.push_back(leptonMatch.gen_particle_lastCopy->charge());
            trainTuple().lepton_genLast_pt.push_back(static_cast<float>(leptonMatch.gen_particle_lastCopy->polarP4().pt()));
            trainTuple().lepton_genLast_eta.push_back(static_cast<float>(leptonMatch.gen_particle_lastCopy->polarP4().eta()));
            trainTuple().lepton_genLast_phi.push_back(static_cast<float>(leptonMatch.gen_particle_lastCopy->polarP4().phi()));
            trainTuple().lepton_genLast_mass.push_back(static_cast<float>(leptonMatch.gen_particle_lastCopy->polarP4().mass()));
            trainTuple().lepton_genFirst_charge.push_back(leptonMatch.gen_particle_firstCopy->charge());
            trainTuple().lepton_genFirst_pt.push_back(static_cast<float>(leptonMatch.gen_particle_firstCopy->polarP4().pt()));
            trainTuple().lepton_genFirst_eta.push_back(static_cast<float>(leptonMatch.gen_particle_firstCopy->polarP4().eta()));
            trainTuple().lepton_genFirst_phi.push_back(static_cast<float>(leptonMatch.gen_particle_firstCopy->polarP4().phi()));
            trainTuple().lepton_genFirst_mass.push_back(static_cast<float>(leptonMatch.gen_particle_firstCopy->polarP4().mass()));

            for(auto daughter : leptonMatch.visible_daughters) {
                trainTuple().lepton_gen_vis_daug_index.push_back(n);
                trainTuple().lepton_gen_vis_daug_pdg.push_back(daughter->pdgId());
                trainTuple().lepton_gen_vis_daug_pt.push_back(static_cast<float>(daughter->polarP4().pt()));
                trainTuple().lepton_gen_vis_daug_eta.push_back(static_cast<float>(daughter->polarP4().eta()));
                trainTuple().lepton_gen_vis_daug_phi.push_back(static_cast<float>(daughter->polarP4().phi()));
                trainTuple().lepton_gen_vis_daug_mass.push_back(static_cast<float>(daughter->polarP4().mass()));
            }
            for(auto daughter : leptonMatch.visible_rad) {
                trainTuple().lepton_gen_vis_rad_index.push_back(n);
                trainTuple().lepton_gen_vis_rad_pdg.push_back(daughter->pdgId());
                trainTuple().lepton_gen_vis_rad_pt.push_back(static_cast<float>(daughter->polarP4().pt()));
                trainTuple().lepton_gen_vis_rad_eta.push_back(static_cast<float>(daughter->polarP4().eta()));
                trainTuple().lepton_gen_vis_rad_phi.push_back(static_cast<float>(daughter->polarP4().phi()));
                trainTuple().lepton_gen_vis_rad_mass.push_back(static_cast<float>(daughter->polarP4().mass()));
            }

            trainTuple().lepton_gen_visible_p4_pt.push_back(static_cast<float>(leptonMatch.visible_p4.pt()));
            trainTuple().lepton_gen_visible_p4_eta.push_back(static_cast<float>(leptonMatch.visible_p4.eta()));
            trainTuple().lepton_gen_visible_p4_phi.push_back(static_cast<float>(leptonMatch.visible_p4.phi()));
            trainTuple().lepton_gen_visible_p4_mass.push_back(static_cast<float>(leptonMatch.visible_p4.mass()));
            trainTuple().lepton_gen_visible_rad_p4_pt.push_back(static_cast<float>(leptonMatch.visible_rad_p4.pt()));
            trainTuple().lepton_gen_visible_rad_p4_eta.push_back(static_cast<float>(leptonMatch.visible_rad_p4.eta()));
            trainTuple().lepton_gen_visible_rad_p4_phi.push_back(static_cast<float>(leptonMatch.visible_rad_p4.phi()));
            trainTuple().lepton_gen_visible_rad_p4_mass.push_back(static_cast<float>(leptonMatch.visible_rad_p4.mass()));
            trainTuple().lepton_gen_n_charged_hadrons.push_back(leptonMatch.n_charged_hadrons);
            trainTuple().lepton_gen_n_neutral_hadrons.push_back(leptonMatch.n_neutral_hadrons);
            trainTuple().lepton_gen_n_gammas.push_back(leptonMatch.n_gammas);
            trainTuple().lepton_gen_n_gammas_rad.push_back(leptonMatch.n_gammas_rad);
        }


    }

    void FillL1Objects(const l1t::TauBxCollection& l1Taus)
    {
        for(auto iter = l1Taus.begin(0); iter != l1Taus.end(0); ++iter) {
            trainTuple().l1Tau_pt.push_back(static_cast<float>(iter->polarP4().pt()));
            trainTuple().l1Tau_eta.push_back(static_cast<float>(iter->polarP4().eta()));
            trainTuple().l1Tau_phi.push_back(static_cast<float>(iter->polarP4().phi()));
            trainTuple().l1Tau_mass.push_back(static_cast<float>(iter->polarP4().mass()));
            trainTuple().l1Tau_hwIso.push_back(iter->hwIso());
            trainTuple().l1Tau_hwQual.push_back(iter->hwQual());
        }
    }

    void FillPixelTracks(const reco::TrackCollection& pixelTracks)
    {
        for(unsigned n = 0; n < pixelTracks.size(); ++n){
            trainTuple().track_pt.push_back(pixelTracks.at(n).pt());
            trainTuple().track_eta.push_back(pixelTracks.at(n).eta());
            trainTuple().track_phi.push_back(pixelTracks.at(n).phi());
            trainTuple().track_outerOk.push_back(pixelTracks.at(n).outerOk());
            trainTuple().track_innerOk.push_back(pixelTracks.at(n).innerOk());
            trainTuple().track_found.push_back(pixelTracks.at(n).found());
            trainTuple().track_lost.push_back(pixelTracks.at(n).lost());
            trainTuple().track_chi2.push_back(pixelTracks.at(n).chi2());
            trainTuple().track_ndof.push_back(pixelTracks.at(n).ndof());
            trainTuple().track_charge.push_back(pixelTracks.at(n).charge());
            trainTuple().track_algo.push_back(static_cast<int>(pixelTracks.at(n).algo()));
            trainTuple().track_qualityMask.push_back(static_cast<unsigned>(pixelTracks.at(n).qualityMask()));
            trainTuple().track_dxy.push_back(pixelTracks.at(n).dxy());
            trainTuple().track_dz.push_back(pixelTracks.at(n).dz());
            trainTuple().track_vx.push_back(pixelTracks.at(n).vx());
            trainTuple().track_vy.push_back(pixelTracks.at(n).vy());
            trainTuple().track_vz.push_back(pixelTracks.at(n).vz());
            trainTuple().track_ptError.push_back(pixelTracks.at(n).ptError());
            trainTuple().track_etaError.push_back(pixelTracks.at(n).etaError());
            trainTuple().track_phiError.push_back(pixelTracks.at(n).phiError());
            trainTuple().track_dxyError.push_back(pixelTracks.at(n).dxyError());
            trainTuple().track_dzError.push_back(pixelTracks.at(n).dzError());
        }
    }

    void FillCaloTowers(const CaloTowerCollection& caloTowers)
    {
        for(auto iter = caloTowers.begin(); iter != caloTowers.end(); ++iter){
            trainTuple().caloTower_pt.push_back(iter->polarP4().pt());
            trainTuple().caloTower_eta.push_back(iter->polarP4().eta());
            trainTuple().caloTower_phi.push_back(iter->polarP4().phi());
            trainTuple().caloTower_energy.push_back(iter->polarP4().energy());
            trainTuple().caloTower_emEnergy.push_back(iter->emEnergy());
            trainTuple().caloTower_hadEnergy.push_back(iter->hadEnergy());
            trainTuple().caloTower_outerEnergy.push_back(iter->outerEnergy());
            trainTuple().caloTower_emPosition_x.push_back(iter->emPosition().x());
            trainTuple().caloTower_emPosition_y.push_back(iter->emPosition().y());
            trainTuple().caloTower_emPosition_z.push_back(iter->emPosition().z());
            trainTuple().caloTower_hadPosition_x.push_back(iter->hadPosition().x());
            trainTuple().caloTower_hadPosition_y.push_back(iter->hadPosition().y());
            trainTuple().caloTower_hadPosition_z.push_back(iter->hadPosition().z());
            trainTuple().caloTower_hadEnergyHeOuterLayer.push_back(iter->hadEnergyHeOuterLayer());
            trainTuple().caloTower_hadEnergyHeInnerLayer.push_back(iter->hadEnergyHeInnerLayer());
            trainTuple().caloTower_energyInHB.push_back(iter->energyInHB());
            trainTuple().caloTower_energyInHE.push_back(iter->energyInHE());
            trainTuple().caloTower_energyInHF.push_back(iter->energyInHF());
            trainTuple().caloTower_energyInHO.push_back(iter->energyInHO());
            trainTuple().caloTower_numBadEcalCells.push_back(iter->numBadEcalCells());
            trainTuple().caloTower_numRecoveredEcalCells.push_back(iter->numRecoveredEcalCells());
            trainTuple().caloTower_numProblematicEcalCells.push_back(iter->numProblematicEcalCells());
            trainTuple().caloTower_numBadHcalCells.push_back(iter->numBadHcalCells());
            trainTuple().caloTower_numRecoveredHcalCells.push_back(iter->numRecoveredHcalCells());
            trainTuple().caloTower_numProblematicHcalCells.push_back(iter->numProblematicHcalCells());
            trainTuple().caloTower_ecalTime.push_back(iter->ecalTime());
            trainTuple().caloTower_hcalTime.push_back(iter->hcalTime());
        }
    }

    void FillCaloTaus(const reco::CaloJetCollection& caloTaus)
    {
        for(const auto& caloTau : caloTaus){
            trainTuple().caloTau_pt.push_back(caloTau.p4().pt());
            trainTuple().caloTau_eta.push_back(caloTau.p4().eta());
            trainTuple().caloTau_phi.push_back(caloTau.p4().phi());
            trainTuple().caloTau_energy.push_back(caloTau.p4().energy());
            trainTuple().caloTau_maxEInEmTowers.push_back(caloTau.maxEInEmTowers());
            trainTuple().caloTau_maxEInHadTowers.push_back(caloTau.maxEInHadTowers());
            trainTuple().caloTau_energyFractionHadronic.push_back(caloTau.energyFractionHadronic());
            trainTuple().caloTau_emEnergyFraction.push_back(caloTau.emEnergyFraction());
            trainTuple().caloTau_hadEnergyInHB.push_back(caloTau.hadEnergyInHB());
            trainTuple().caloTau_hadEnergyInHO.push_back(caloTau.hadEnergyInHO());
            trainTuple().caloTau_hadEnergyInHE.push_back(caloTau.hadEnergyInHE());
            trainTuple().caloTau_hadEnergyInHF.push_back(caloTau.hadEnergyInHF());
            trainTuple().caloTau_emEnergyInEB.push_back(caloTau.emEnergyInEB());
            trainTuple().caloTau_emEnergyInEE.push_back(caloTau.emEnergyInEE());
            trainTuple().caloTau_emEnergyInHF.push_back(caloTau.emEnergyInHF());
            trainTuple().caloTau_towersArea.push_back(caloTau.towersArea());
            trainTuple().caloTau_n90.push_back(caloTau.n90());
            trainTuple().caloTau_n60.push_back(caloTau.n60());
        }
    }

    void FillPFCandidates(const std::vector<reco::PFCandidate>& cands)
    {
        for(const auto& cand : cands) {
            trainTuple().pfCand_pt.push_back(static_cast<float>(cand.polarP4().pt()));
            trainTuple().pfCand_eta.push_back(static_cast<float>(cand.polarP4().eta()));
            trainTuple().pfCand_phi.push_back(static_cast<float>(cand.polarP4().phi()));
            trainTuple().pfCand_mass.push_back(static_cast<float>(cand.polarP4().mass()));

            trainTuple().pfCand_pvAssociationQuality.push_back(default_int_value);
            trainTuple().pfCand_fromPV.push_back(default_int_value);
            trainTuple().pfCand_puppiWeight.push_back(default_value);
            trainTuple().pfCand_puppiWeightNoLep.push_back(default_value);
            trainTuple().pfCand_pdgId.push_back(cand.pdgId());
            trainTuple().pfCand_charge.push_back(cand.charge());
            int lostInnerHits = cand.bestTrack() != nullptr ? cand.bestTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) : default_int_value;
            trainTuple().pfCand_lostInnerHits.push_back(lostInnerHits);
            int numberOfPixelHits = cand.bestTrack() != nullptr ? cand.bestTrack()->hitPattern().numberOfValidPixelHits() : default_int_value;
            trainTuple().pfCand_numberOfPixelHits.push_back(numberOfPixelHits);

            trainTuple().pfCand_vertex_x.push_back(static_cast<float>(cand.vertex().x()));
            trainTuple().pfCand_vertex_y.push_back(static_cast<float>(cand.vertex().y()));
            trainTuple().pfCand_vertex_z.push_back(static_cast<float>(cand.vertex().z()));

            const bool hasTrackDetails = cand.bestTrack() != nullptr;
            trainTuple().pfCand_hasTrackDetails.push_back(hasTrackDetails);
            trainTuple().pfCand_dxy.push_back(hasTrackDetails ? cand.bestTrack()->dxy() : default_value);
            trainTuple().pfCand_dxy_error.push_back(hasTrackDetails ? cand.dxyError() : default_value);
            trainTuple().pfCand_dz.push_back(hasTrackDetails ? cand.bestTrack()->dz() : default_value);
            trainTuple().pfCand_dz_error.push_back(hasTrackDetails ? cand.dzError() : default_value);
            trainTuple().pfCand_track_chi2.push_back(
                        hasTrackDetails ? static_cast<float>(cand.bestTrack()->chi2()) : default_value);
            trainTuple().pfCand_track_ndof.push_back(
                        hasTrackDetails ? static_cast<float>(cand.bestTrack()->ndof()) : default_value);

            float hcal_fraction = cand.rawHcalEnergy()/(cand.rawHcalEnergy()+cand.rawEcalEnergy());
            trainTuple().pfCand_hcalFraction.push_back(hcal_fraction);
            trainTuple().pfCand_rawCaloFraction.push_back((cand.rawEcalEnergy()+cand.rawHcalEnergy())/cand.energy());

            //electron vars
            trainTuple().pfCand_ele_trackMomentumAtVtx.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->trackMomentumAtVtx().R() : default_value);
            trainTuple().pfCand_ele_trackMomentumAtCalo.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->trackMomentumAtCalo().R() : default_value);
            trainTuple().pfCand_ele_trackMomentumOut.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->trackMomentumOut().R() : default_value);
            trainTuple().pfCand_ele_trackMomentumAtEleClus.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->trackMomentumAtEleClus().R() : default_value);
            trainTuple().pfCand_ele_trackMomentumAtVtxWithConstraint.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->trackMomentumAtVtxWithConstraint().R() : default_value);
            trainTuple().pfCand_ele_ecalEnergy.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->ecalEnergy() : default_value);
            trainTuple().pfCand_ele_ecalEnergy_error.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->ecalEnergyError() : default_value);
            trainTuple().pfCand_ele_eSuperClusterOverP.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->eSuperClusterOverP() : default_value);
            trainTuple().pfCand_ele_eSeedClusterOverP.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->eSeedClusterOverP() : default_value);
            trainTuple().pfCand_ele_eSeedClusterOverPout.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->eSeedClusterOverPout() : default_value);
            trainTuple().pfCand_ele_eEleClusterOverPout.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->eEleClusterOverPout() : default_value);
            trainTuple().pfCand_ele_deltaEtaSuperClusterTrackAtVtx.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaEtaSuperClusterTrackAtVtx() : default_value);
            trainTuple().pfCand_ele_deltaEtaSeedClusterTrackAtCalo.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaEtaSeedClusterTrackAtCalo() : default_value);
            trainTuple().pfCand_ele_deltaEtaEleClusterTrackAtCalo.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaEtaEleClusterTrackAtCalo() : default_value);
            trainTuple().pfCand_ele_deltaPhiEleClusterTrackAtCalo.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaPhiEleClusterTrackAtCalo() : default_value);
            trainTuple().pfCand_ele_deltaPhiSuperClusterTrackAtVtx.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaPhiSuperClusterTrackAtVtx() : default_value);
            trainTuple().pfCand_ele_deltaPhiSeedClusterTrackAtCalo.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->deltaPhiSeedClusterTrackAtCalo() : default_value);

            trainTuple().pfCand_ele_mvaInput_earlyBrem.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->mvaInput().earlyBrem : default_value);
            trainTuple().pfCand_ele_mvaInput_lateBrem.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->mvaInput().lateBrem : default_value);
            trainTuple().pfCand_ele_mvaInput_sigmaEtaEta.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->mvaInput().sigmaEtaEta : default_value);
            trainTuple().pfCand_ele_mvaInput_hadEnergy.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->mvaInput().hadEnergy : default_value);
            trainTuple().pfCand_ele_mvaInput_deltaEta.push_back(cand.gsfElectronRef().isNonnull() ? cand.gsfElectronRef()->mvaInput().deltaEta : default_value);

            bool hasClosest = cand.gsfElectronRef().isNonnull() && cand.gsfElectronRef()->closestCtfTrackRef().isNonnull();
            trainTuple().pfCand_ele_closestCtfTrack_normalizedChi2.push_back(hasClosest
                        ? static_cast<float>(cand.gsfElectronRef()->closestCtfTrackRef()->normalizedChi2()) : default_value);
            trainTuple().pfCand_ele_closestCtfTrack_numberOfValidHits.push_back(
                        hasClosest ? cand.gsfElectronRef()->closestCtfTrackRef()->numberOfValidHits() : default_int_value);
        }
    }

    void FillElectrons(const std::vector<reco::RecoEcalCandidate>& electrons)
    {
        for(const auto& ele : electrons) {
            trainTuple().ele_pt.push_back(static_cast<float>(ele.polarP4().pt()));
            trainTuple().ele_eta.push_back(static_cast<float>(ele.polarP4().eta()));
            trainTuple().ele_phi.push_back(static_cast<float>(ele.polarP4().phi()));
            trainTuple().ele_mass.push_back(static_cast<float>(ele.polarP4().mass()));
            float cc_ele_energy, cc_gamma_energy;
            int cc_n_gamma;
            CalculateElectronClusterVars(ele, cc_ele_energy, cc_gamma_energy, cc_n_gamma);
            trainTuple().ele_cc_ele_energy.push_back(cc_ele_energy);
            trainTuple().ele_cc_gamma_energy.push_back(cc_gamma_energy);
            trainTuple().ele_cc_n_gamma.push_back(cc_n_gamma);

            const auto& gsfTrack = ele.gsfTrack();
            trainTuple().ele_gsfTrack_normalizedChi2.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->normalizedChi2()) : default_value);
            trainTuple().ele_gsfTrack_numberOfValidHits.push_back(
                        gsfTrack.isNonnull() ? gsfTrack->numberOfValidHits() : default_int_value);
            trainTuple().ele_gsfTrack_pt.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->pt()) : default_value);
            trainTuple().ele_gsfTrack_pt_error.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->ptError()) : default_value);

        }
    }

    void FillMuons(const std::vector<reco::Muon>& muons, const std::vector<reco::PFCandidate>& cands)
    {
        for(const auto& muon : muons) {
            trainTuple().muon_pt.push_back(static_cast<float>(muon.polarP4().pt()));
            trainTuple().muon_eta.push_back(static_cast<float>(muon.polarP4().eta()));
            trainTuple().muon_phi.push_back(static_cast<float>(muon.polarP4().phi()));
            trainTuple().muon_mass.push_back(static_cast<float>(muon.polarP4().mass()));
            bool hasBestTrack = muon.bestTrack() != nullptr;
            trainTuple().muon_dxy.push_back(hasBestTrack ? static_cast<float>(muon.bestTrack()->dxy()) : default_value);
            trainTuple().muon_dxy_error.push_back(hasBestTrack ? static_cast<float>(std::abs(muon.bestTrack()->dxy())/muon.dxyError()) : default_value);

            const bool normalizedChi2_valid = hasBestTrack && muon.bestTrack()->outerOk() && muon.bestTrack()->normalizedChi2() >= 0;
            if(normalizedChi2_valid){
                trainTuple().muon_normalizedChi2.push_back(muon.bestTrack()->normalizedChi2());
                if(muon.bestTrack()->innerOk())
                    trainTuple().muon_numberOfValidHits.push_back(static_cast<int>(muon.bestTrack()->numberOfValidHits()));
            }

            double segmentCompatibility = muon::segmentCompatibility(muon,reco::Muon::SegmentAndTrackArbitration);
            trainTuple().muon_segmentCompatibility.push_back(static_cast<float>(segmentCompatibility));
            trainTuple().muon_caloCompatibility.push_back(muon.caloCompatibility());

            bool pfEcalEnergy_valid = false;
            double pfEcalEnergy = 0.0;
            for (const reco::PFCandidate& pfcand : cands) {
                if (pfcand.muonRef().isNonnull()) {
                    if(&muon == &(*pfcand.muonRef())){
                        pfEcalEnergy_valid = pfcand.ecalEnergy() >= 0;
                        pfEcalEnergy = pfcand.ecalEnergy();
                    }
                }
            }

            trainTuple().muon_pfEcalEnergy.push_back(pfEcalEnergy_valid ? pfEcalEnergy : default_value);

            const MuonHitMatch hit_match(muon);
            for(int subdet : MuonHitMatch::ConsideredSubdets()) {
                const std::string& subdetName = MuonHitMatch::SubdetName(subdet);
                for(int station = MuonHitMatch::first_station_id; station <= MuonHitMatch::last_station_id; ++station) {
                    const std::string matches_branch_name = "muon_n_matches_" + subdetName + "_"
                            + std::to_string(station);
                    const std::string hits_branch_name = "muon_n_hits_" + subdetName + "_" + std::to_string(station);

                    const unsigned n_matches = hit_match.NMatches(subdet, station);
                    const unsigned n_hits = hit_match.NHits(subdet, station);
                    trainTuple.get<std::vector<int>>(matches_branch_name).push_back(static_cast<int>(n_matches));
                    trainTuple.get<std::vector<int>>(hits_branch_name).push_back(static_cast<int>(n_hits));
                }
            }
        }
    }

    static float CalculateGottfriedJacksonAngleDifference(const reco::PFTau& tau)
    {
        double gj_diff;
        if(::tau_analysis::CalculateGottfriedJacksonAngleDifference(tau, gj_diff))
            return static_cast<float>(gj_diff);
        return default_value;
    }

    static void CalculateElectronClusterVars(const reco::RecoEcalCandidate& ele, float& cc_ele_energy, float& cc_gamma_energy,
                                             int& cc_n_gamma)
    {
        cc_ele_energy = cc_gamma_energy = 0;
        cc_n_gamma = 0;
        const auto& superCluster = ele.superCluster();
        if(superCluster.isNonnull() && superCluster.isAvailable() && superCluster->clusters().isNonnull()
                && superCluster->clusters().isAvailable()) {
            for(auto iter = superCluster->clustersBegin(); iter != superCluster->clustersEnd(); ++iter) {
                const float energy = static_cast<float>((*iter)->energy());
                if(iter == superCluster->clustersBegin())
                    cc_ele_energy += energy;
                else {
                    cc_gamma_energy += energy;
                    ++cc_n_gamma;
                }
            }
        } else {
            cc_ele_energy = cc_gamma_energy = default_value;
            cc_n_gamma = default_int_value;
        }
    }

private:
    const bool isMC, storeJetsWithoutTau, requireGenMatch;
    TauJetBuilderSetup builderSetup;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<double> rho_token;
    edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> electrons_token;
    edm::EDGetTokenT<reco::MuonCollection> muons_token;
    edm::EDGetTokenT<std::vector<reco::PFTau>> taus_token;
    edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;
    edm::EDGetTokenT<CaloTowerCollection> caloTowers_token;
    edm::EDGetTokenT<reco::TrackCollection> pixelTracks_token;
    edm::EDGetTokenT<reco::CaloJetCollection> caloTaus_token;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> cands_token;
    edm::EDGetTokenT<reco::PFTauDiscriminator> decayMode_token;
    edm::EDGetTokenT<TauDiscriminator> chargedIsoPtSum_inputToken;
    edm::EDGetTokenT<TauDiscriminator> chargedIsoPtSumdR03_inputToken;
    edm::EDGetTokenT<TauDiscriminator> neutralIsoPtSum_inputToken;
    edm::EDGetTokenT<TauDiscriminator> neutralIsoPtSumdR03_inputToken;
    edm::EDGetTokenT<TauDiscriminator> footprintCorrection_inputToken;
    edm::EDGetTokenT<TauDiscriminator> footprintCorrectiondR03_inputToken;
    edm::EDGetTokenT<TauDiscriminator> neutralIsoPtSumWeight_inputToken;
    edm::EDGetTokenT<TauDiscriminator> neutralIsoPtSumWeightdR03_inputToken;
    edm::EDGetTokenT<TauDiscriminator> photonPtSumOutsideSignalCone_inputToken;
    edm::EDGetTokenT<TauDiscriminator> photonPtSumOutsideSignalConedR03_inputToken;
    edm::EDGetTokenT<TauDiscriminator> puCorrPtSum_inputToken;
    edm::EDGetTokenT<edm::AssociationVector<reco::PFTauRefProd, std::vector<reco::PFTauTransverseImpactParameterRef>>> PFTauTransverseImpactParameters_token;
    const edm::EDGetTokenT<TauDiscriminator> deepTauVSe_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> deepTauVSmu_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> deepTauVSjet_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> looseIsoAbs_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> looseIsoRel_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> mediumIsoAbs_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> mediumIsoRel_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> tightIsoAbs_inputToken;
    const edm::EDGetTokenT<TauDiscriminator> tightIsoRel_inputToken;

    TrainTupleProducerData* data;
    train_tuple::TrainTuple& trainTuple;
    tau_tuple::SummaryTuple& summaryTuple;
};

} // namespace tau_analysis

#include "FWCore/Framework/interface/MakerMacros.h"
using TrainTupleProducer = tau_analysis::TrainTupleProducer;
DEFINE_FWK_MODULE(TrainTupleProducer);
