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
#include "TauML/Analysis/include/TauTuple.h"
#include "TauML/Analysis/include/SummaryTuple.h"
#include "TauML/Analysis/include/TauIdResults.h"
#include "TauML/Production/include/GenTruthTools.h"
#include "TauML/Production/include/TauAnalysis.h"
#include "TauML/Production/include/MuonHitMatch.h"
#include "TauML/Production/include/TauJet.h"

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

namespace tau_analysis {

struct TauTupleProducerData {
    using clock = std::chrono::system_clock;

    const clock::time_point start;
    tau_tuple::TauTuple tauTuple;
    tau_tuple::SummaryTuple summaryTuple;
    std::mutex mutex;

private:
    size_t n_producers;

    TauTupleProducerData(TFile& file) :
        start(clock::now()),
        tauTuple("taus", &file, false),
        summaryTuple("summary", &file, false),
        n_producers(0)
    {
        summaryTuple().numberOfProcessedEvents = 0;
    }

    ~TauTupleProducerData() {}

public:

    static TauTupleProducerData* RequestGlobalData()
    {
        TauTupleProducerData* data = GetGlobalData();
        if(data == nullptr)
            throw cms::Exception("TauTupleProducerData") << "Request after all data copies were released.";
        {
            std::lock_guard<std::mutex> lock(data->mutex);
            ++data->n_producers;
            std::cout << "New request of TauTupleProducerData. Total number of producers = " << data->n_producers
                      << "." << std::endl;
        }
        return data;
    }

    static void ReleaseGlobalData()
    {
        TauTupleProducerData*& data = GetGlobalData();
        if(data == nullptr)
            throw cms::Exception("TauTupleProducerData") << "Another release after all data copies were released.";
        {
            std::lock_guard<std::mutex> lock(data->mutex);
            if(!data->n_producers)
                throw cms::Exception("TauTupleProducerData") << "Release before any request.";
            --data->n_producers;
            std::cout << "TauTupleProducerData has been released. Total number of producers = " << data->n_producers
                      << "." << std::endl;
            if(!data->n_producers) {
                data->tauTuple.Write();
                const auto stop = clock::now();
                data->summaryTuple().exeTime = static_cast<unsigned>(
                            std::chrono::duration_cast<std::chrono::seconds>(stop - data->start).count());
                data->summaryTuple.Fill();
                data->summaryTuple.Write();
                delete data;
                data = nullptr;
                std::cout << "TauTupleProducerData has been destroyed." << std::endl;
            }
        }

    }

private:
    static TauTupleProducerData*& GetGlobalData()
    {
        static TauTupleProducerData* data = InitializeGlobalData();
        return data;
    }

    static TauTupleProducerData* InitializeGlobalData()
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kZLIB);
        file.SetCompressionLevel(9);
        TauTupleProducerData* data = new TauTupleProducerData(file);
        std::cout << "TauTupleProducerData has been created." << std::endl;
        return data;
    }
};

class TauTupleProducer : public edm::EDAnalyzer {
public:
    using TauDiscriminator = reco::PFTauDiscriminator;

    TauTupleProducer(const edm::ParameterSet& cfg) :
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
        //jets_token(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
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
        data(TauTupleProducerData::RequestGlobalData()),
        tauTuple(data->tauTuple),
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
    static constexpr float default_value = tau_tuple::DefaultFillValue<float>();
    static constexpr int default_int_value = tau_tuple::DefaultFillValue<int>();

    virtual void analyze(const edm::Event& event, const edm::EventSetup&) override
    {
        std::lock_guard<std::mutex> lock(data->mutex);
        summaryTuple().numberOfProcessedEvents++;

        tauTuple().run  = event.id().run();
        tauTuple().lumi = event.id().luminosityBlock();
        tauTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::Vertex>> vertices;
        event.getByToken(vertices_token, vertices);
        tauTuple().npv = static_cast<int>(vertices->size());
        edm::Handle<double> rho;
        event.getByToken(rho_token, rho);
        tauTuple().rho = static_cast<float>(*rho);

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            tauTuple().genEventWeight = static_cast<float>(genEvent->weight());

            edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
            event.getByToken(puInfo_token, puInfo);
            tauTuple().npu = gen_truth::GetNumberOfPileUpInteractions(puInfo);
        }

        const auto& PV = vertices->at(0);
        tauTuple().pv_x = static_cast<float>(PV.position().x());
        tauTuple().pv_y = static_cast<float>(PV.position().y());
        tauTuple().pv_z = static_cast<float>(PV.position().z());
        tauTuple().pv_chi2 = static_cast<float>(PV.chi2());
        tauTuple().pv_ndof = static_cast<float>(PV.ndof());

        edm::Handle<std::vector<reco::RecoEcalCandidate>> electrons;
        event.getByToken(electrons_token, electrons);

        edm::Handle<reco::MuonCollection> muons;
        event.getByToken(muons_token, muons);

        edm::Handle<std::vector<reco::PFTau>> taus;
        event.getByToken(taus_token, taus);

        // edm::Handle<pat::JetCollection> jets;
        // event.getByToken(jets_token, jets);
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

        TauJetBuilder builder(builderSetup, jets, *taus, *cands, *electrons, *muons, genParticles);
        const auto tauJets = builder.Build();

        for(const TauJet& tauJet : tauJets) {
            const bool has_jet = tauJet.jetIndex >= 0;
            const bool has_tau = tauJet.tauIndex >= 0;
            if(!has_tau && !storeJetsWithoutTau) continue;

            const auto& leptonGenMatch = has_tau ? tauJet.tauGenLeptonMatchResult : tauJet.jetGenLeptonMatchResult;
            const auto& qcdGenMatch = has_tau ? tauJet.tauGenQcdMatchResult : tauJet.jetGenQcdMatchResult;

            if(requireGenMatch && leptonGenMatch.match == GenLeptonMatch::NoMatch
                    && qcdGenMatch.match == GenQcdMatch::NoMatch) continue;

            tauTuple().jet_index = tauJet.jetIndex;
            tauTuple().jet_pt = has_jet ? static_cast<float>(tauJet.jet->p4().pt()) : default_value;
            tauTuple().jet_eta = has_jet ? static_cast<float>(tauJet.jet->p4().eta()) : default_value;
            tauTuple().jet_phi = has_jet ? static_cast<float>(tauJet.jet->p4().phi()) : default_value;
            tauTuple().jet_mass = has_jet ? static_cast<float>(tauJet.jet->p4().mass()) : default_value;
            boost::optional<pat::Jet> uncorrected_jet;
            if(has_jet)
                uncorrected_jet = tauJet.jet->correctedJet("Uncorrected");
            tauTuple().jet_neutralHadronEnergyFraction = has_jet
                    ? uncorrected_jet->neutralHadronEnergyFraction() : default_value;
            tauTuple().jet_neutralEmEnergyFraction = has_jet
                    ? uncorrected_jet->neutralEmEnergyFraction() : default_value;
            tauTuple().jet_nConstituents = has_jet ? uncorrected_jet->nConstituents() : default_int_value;
            tauTuple().jet_chargedMultiplicity = has_jet ? uncorrected_jet->chargedMultiplicity() : default_int_value;
            tauTuple().jet_neutralMultiplicity = has_jet ? uncorrected_jet->neutralMultiplicity() : default_int_value;
            tauTuple().jet_partonFlavour = has_jet ? tauJet.jet->partonFlavour() : default_int_value;
            tauTuple().jet_hadronFlavour = has_jet ? tauJet.jet->hadronFlavour() : default_int_value;
            const reco::GenJet* genJet = has_jet ? tauJet.jet->genJet() : nullptr;
            tauTuple().jet_has_gen_match = genJet != nullptr;
            tauTuple().jet_gen_pt = genJet != nullptr ? static_cast<float>(genJet->polarP4().pt()) : default_value;
            tauTuple().jet_gen_eta = genJet != nullptr ? static_cast<float>(genJet->polarP4().eta()) : default_value;
            tauTuple().jet_gen_phi = genJet != nullptr ? static_cast<float>(genJet->polarP4().phi()) : default_value;
            tauTuple().jet_gen_mass = genJet != nullptr ? static_cast<float>(genJet->polarP4().mass()) : default_value;
            tauTuple().jet_gen_n_b = genJet != nullptr
                    ? static_cast<int>(tauJet.jet->jetFlavourInfo().getbHadrons().size()) : default_int_value;
            tauTuple().jet_gen_n_c = genJet != nullptr
                    ? static_cast<int>(tauJet.jet->jetFlavourInfo().getcHadrons().size()) : default_int_value;


            const reco::PFTau* tau = tauJet.tau;

            tauTuple().jetTauMatch = static_cast<int>(tauJet.jetTauMatch);
            tauTuple().tau_index = tauJet.tauIndex;
            tauTuple().tau_pt = has_tau ? static_cast<float>(tau->polarP4().pt()) : default_value;
            tauTuple().tau_eta = has_tau ? static_cast<float>(tau->polarP4().eta()) : default_value;
            tauTuple().tau_phi = has_tau ? static_cast<float>(tau->polarP4().phi()) : default_value;
            tauTuple().tau_mass = has_tau ? static_cast<float>(tau->polarP4().mass()) : default_value;
            tauTuple().tau_charge = has_tau ? tau->charge() : default_int_value;

            FillGenMatchResult(leptonGenMatch, qcdGenMatch);

            tauTuple().tau_decayMode = has_tau ? tau->decayMode() : default_int_value;
            // tauTuple().tau_decayModeFinding = has_tau ? tau->tauID("decayModeFinding") > 0.5f : default_int_value;
            tauTuple().tau_decayModeFindingNewDMs = has_tau ? decayModesNew->value(tauJet.tauIndex)
                                                            : default_int_value;
            tauTuple().tau_chargedIsoPtSum = has_tau ? chargedIsoPtSum->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_chargedIsoPtSumdR03 = has_tau ? chargedIsoPtSumdR03->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_footprintCorrection = has_tau ? footprintCorrection->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_footprintCorrectiondR03 = has_tau ? footprintCorrectiondR03->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_neutralIsoPtSum = has_tau ? neutralIsoPtSum->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_neutralIsoPtSumWeight = has_tau ? neutralIsoPtSumWeight->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_neutralIsoPtSumWeightdR03 = has_tau ? neutralIsoPtSumWeightdR03->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_neutralIsoPtSumdR03 = has_tau ? neutralIsoPtSumdR03->value(tauJet.tauIndex) : default_value;
            tauTuple().tau_photonPtSumOutsideSignalCone = has_tau ? photonPtSumOutsideSignalCone->value(tauJet.tauIndex)
                                                              : default_value;
            tauTuple().tau_photonPtSumOutsideSignalConedR03 = has_tau ? photonPtSumOutsideSignalConedR03->value(tauJet.tauIndex)
                                                                  : default_value;
            tauTuple().tau_puCorrPtSum = has_tau ? puCorrPtSum->value(tauJet.tauIndex) : default_value;

            analysis::TauIdResults cutBasedRelIso;
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Loose,looseIsoRel->value(tauJet.tauIndex) > 0.5);
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Medium,mediumIsoRel->value(tauJet.tauIndex) > 0.5);
            cutBasedRelIso.SetResult(analysis::DiscriminatorWP::Tight,tightIsoRel->value(tauJet.tauIndex) > 0.5);
            tauTuple().tau_cutBasedRelIso = cutBasedRelIso.GetResultBits();

            analysis::TauIdResults cutBasedAbsIso;
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Loose,looseIsoAbs->value(tauJet.tauIndex) > 0.5);
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Medium,mediumIsoAbs->value(tauJet.tauIndex) > 0.5);
            cutBasedAbsIso.SetResult(analysis::DiscriminatorWP::Tight,tightIsoAbs->value(tauJet.tauIndex) > 0.5);
            tauTuple().tau_cutBasedAbsIso = cutBasedAbsIso.GetResultBits();

            tauTuple().tau_byDeepTau2017v2VSeraw = static_cast<float>(deepTau_VSe->value(tauJet.tauIndex));
            tauTuple().tau_byDeepTau2017v2VSmuraw = static_cast<float>(deepTau_VSmu->value(tauJet.tauIndex));
            tauTuple().tau_byDeepTau2017v2VSjetraw = static_cast<float>(deepTau_VSjet->value(tauJet.tauIndex));

            // togli da ntuple
            // tauTuple().tau_dxy_pca_x = has_tau ? tau->dxy_PCA().x() : default_value;
            // tauTuple().tau_dxy_pca_y = has_tau ? tau->dxy_PCA().y() : default_value;
            // tauTuple().tau_dxy_pca_z = has_tau ? tau->dxy_PCA().z() : default_value;

            auto impactParam = PFTauTransverseImpactParameters->value(tauJet.tauIndex);

            tauTuple().tau_dxy = has_tau ? impactParam->dxy() : default_value;
            tauTuple().tau_dxy_error = has_tau ? impactParam->dxy_error() : default_value;
            tauTuple().tau_ip3d = has_tau ? impactParam->ip3d() : default_value;
            tauTuple().tau_ip3d_error = has_tau ? impactParam->ip3d_error() : default_value;
            const bool has_sv = has_tau && impactParam->hasSecondaryVertex();
            tauTuple().tau_hasSecondaryVertex = has_tau ? impactParam->hasSecondaryVertex() : default_int_value;
            tauTuple().tau_sv_x = has_sv ? impactParam->secondaryVertexPos().x() : default_value;
            tauTuple().tau_sv_y = has_sv ? impactParam->secondaryVertexPos().y() : default_value;
            tauTuple().tau_sv_z = has_sv ? impactParam->secondaryVertexPos().z() : default_value;
            tauTuple().tau_flightLength_x = has_tau ? impactParam->flightLength().x() : default_value;
            tauTuple().tau_flightLength_y = has_tau ? impactParam->flightLength().y() : default_value;
            tauTuple().tau_flightLength_z = has_tau ? impactParam->flightLength().z() : default_value;
            tauTuple().tau_flightLength_sig = has_tau ? impactParam->flightLengthSig() : default_value;

            const reco::PFCandidate* leadChargedHadrCand =
                    has_tau ? dynamic_cast<const reco::PFCandidate*>(tau->leadChargedHadrCand().get()) : nullptr;
            tauTuple().tau_dz = leadChargedHadrCand  && leadChargedHadrCand->bestTrack() != nullptr ? leadChargedHadrCand->bestTrack()->dz() : default_value;
            tauTuple().tau_dz_error = leadChargedHadrCand ? leadChargedHadrCand->dzError() : default_value;

            tauTuple().tau_pt_weighted_deta_strip =
                    has_tau ? reco::tau::pt_weighted_deta_strip(*tau, tau->decayMode()) : default_value;
            tauTuple().tau_pt_weighted_dphi_strip =
                    has_tau ? reco::tau::pt_weighted_dphi_strip(*tau, tau->decayMode()) : default_value;
            tauTuple().tau_pt_weighted_dr_signal =
                    has_tau ? reco::tau::pt_weighted_dr_signal(*tau, tau->decayMode()) : default_value;
            tauTuple().tau_pt_weighted_dr_iso =
                    has_tau ? reco::tau::pt_weighted_dr_iso(*tau, tau->decayMode()) : default_value;
            tauTuple().tau_leadingTrackNormChi2 = has_tau ? reco::tau::lead_track_chi2(*tau) : default_value;
            tauTuple().tau_e_ratio = has_tau ? reco::tau::eratio(*tau) : default_value;
            tauTuple().tau_gj_angle_diff = has_tau ? CalculateGottfriedJacksonAngleDifference(*tau) : default_value;
            tauTuple().tau_n_photons =
                    has_tau ? static_cast<int>(reco::tau::n_photons_total(*tau)) : default_int_value;

            float emFraction = -1.;
            float myHCALenergy = 0.;
            float myECALenergy = 0.;
            if(leadChargedHadrCand && leadChargedHadrCand->bestTrack() != nullptr){
                for (const auto& isoPFCand : tau->isolationPFCands()) {
                      myHCALenergy += isoPFCand->hcalEnergy();
                      myECALenergy += isoPFCand->ecalEnergy();
                }
                for (const auto& signalPFCand : tau->signalPFCands()) {
                    myHCALenergy += signalPFCand->hcalEnergy();
                    myECALenergy += signalPFCand->ecalEnergy();
                }
                if (myHCALenergy + myECALenergy != 0.) {
                  emFraction = myECALenergy / (myHCALenergy + myECALenergy);
                }
            }
            tauTuple().tau_emFraction = has_tau ? emFraction : default_value;
            tauTuple().tau_inside_ecal_crack = has_tau ? IsInEcalCrack(tau->p4().Eta()) : default_value;

            const std::vector<reco::CandidatePtr>& signalCands = tau->signalCands();
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
            tauTuple().leadChargedCand_etaAtEcalEntrance =
                    has_tau ? leadChargedCandEtaAtEcalEntrance : default_value;

            FillPFCandidates(tauJet.cands);
            FillElectrons(tauJet.electrons);
            FillMuons(tauJet.muons,tauJet.cands);

            tauTuple.Fill();
        }
    }

    virtual void endJob() override
    {
        TauTupleProducerData::ReleaseGlobalData();
    }

private:
    // static bool PrintTauIdNames(const reco::PFTau& tau)
    // {
    //     static const std::string header(40, '-');
    //
    //     std::set<std::string> tauId_names;
    //     for(const auto& id : tau.tauIDs())
    //         tauId_names.insert(id.first);
    //     std::cout << "Tau IDs:\n" << header << "\n";
    //     for(const std::string& name : tauId_names)
    //         std::cout << name << "\n";
    //     std::cout << header << std::endl;
    //
    //     return true;
    // }

    void FillGenMatchResult(const gen_truth::LeptonMatchResult& leptonMatch, const gen_truth::QcdMatchResult& qcdMatch)
    {
        const bool has_lepton = leptonMatch.match != GenLeptonMatch::NoMatch;
        tauTuple().lepton_gen_match = static_cast<int>(leptonMatch.match);
        tauTuple().lepton_gen_charge = has_lepton ? leptonMatch.gen_particle->charge() : default_int_value;
        tauTuple().lepton_gen_pt = has_lepton ? static_cast<float>(leptonMatch.gen_particle->polarP4().pt())
                                              : default_value;
        tauTuple().lepton_gen_eta = has_lepton ? static_cast<float>(leptonMatch.gen_particle->polarP4().eta())
                                               : default_value;
        tauTuple().lepton_gen_phi = has_lepton ? static_cast<float>(leptonMatch.gen_particle->polarP4().phi())
                                               : default_value;
        tauTuple().lepton_gen_mass = has_lepton ? static_cast<float>(leptonMatch.gen_particle->polarP4().mass())
                                                : default_value;
        for(auto daughter : leptonMatch.visible_daughters) {
            tauTuple().lepton_gen_vis_pdg.push_back(daughter->pdgId());
            tauTuple().lepton_gen_vis_pt.push_back(static_cast<float>(daughter->polarP4().pt()));
            tauTuple().lepton_gen_vis_eta.push_back(static_cast<float>(daughter->polarP4().eta()));
            tauTuple().lepton_gen_vis_phi.push_back(static_cast<float>(daughter->polarP4().phi()));
            tauTuple().lepton_gen_vis_mass.push_back(static_cast<float>(daughter->polarP4().mass()));
        }

        const bool has_qcd = qcdMatch.match != GenQcdMatch::NoMatch;
        tauTuple().qcd_gen_match = static_cast<int>(qcdMatch.match);
        tauTuple().qcd_gen_charge = has_qcd ? qcdMatch.gen_particle->charge() : default_int_value;
        tauTuple().qcd_gen_pt = has_qcd ? static_cast<float>(qcdMatch.gen_particle->polarP4().pt()) : default_value;
        tauTuple().qcd_gen_eta = has_qcd ? static_cast<float>(qcdMatch.gen_particle->polarP4().eta()) : default_value;
        tauTuple().qcd_gen_phi = has_qcd ? static_cast<float>(qcdMatch.gen_particle->polarP4().phi()) : default_value;
        tauTuple().qcd_gen_mass = has_qcd ? static_cast<float>(qcdMatch.gen_particle->polarP4().mass()) : default_value;
    }

    void FillPFCandidates(const std::vector<PFCandDesc>& cands)
    {
        for(const PFCandDesc& cand_desc : cands) {
            const reco::PFCandidate* cand = cand_desc.candidate;

            tauTuple().pfCand_jetDaughter.push_back(cand_desc.jetDaughter);
            tauTuple().pfCand_tauSignal.push_back(cand_desc.tauSignal);
            tauTuple().pfCand_leadChargedHadrCand.push_back(cand_desc.leadChargedHadrCand);
            tauTuple().pfCand_tauIso.push_back(cand_desc.tauIso);

            tauTuple().pfCand_pt.push_back(static_cast<float>(cand->polarP4().pt()));
            tauTuple().pfCand_eta.push_back(static_cast<float>(cand->polarP4().eta()));
            tauTuple().pfCand_phi.push_back(static_cast<float>(cand->polarP4().phi()));
            tauTuple().pfCand_mass.push_back(static_cast<float>(cand->polarP4().mass()));

            tauTuple().pfCand_pvAssociationQuality.push_back(default_int_value);
            tauTuple().pfCand_fromPV.push_back(default_int_value);
            tauTuple().pfCand_puppiWeight.push_back(default_value);
            tauTuple().pfCand_puppiWeightNoLep.push_back(default_value);
            tauTuple().pfCand_pdgId.push_back(cand->pdgId());
            tauTuple().pfCand_charge.push_back(cand->charge());
            int lostInnerHits = cand->bestTrack() != nullptr ? cand->bestTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) : default_int_value;
            tauTuple().pfCand_lostInnerHits.push_back(lostInnerHits);
            int numberOfPixelHits = cand->bestTrack() != nullptr ? cand->bestTrack()->hitPattern().numberOfValidPixelHits() : default_int_value;
            tauTuple().pfCand_numberOfPixelHits.push_back(numberOfPixelHits);

            tauTuple().pfCand_vertex_x.push_back(static_cast<float>(cand->vertex().x()));
            tauTuple().pfCand_vertex_y.push_back(static_cast<float>(cand->vertex().y()));
            tauTuple().pfCand_vertex_z.push_back(static_cast<float>(cand->vertex().z()));

            const bool hasTrackDetails = cand->bestTrack() != nullptr;
            tauTuple().pfCand_hasTrackDetails.push_back(hasTrackDetails);
            tauTuple().pfCand_dxy.push_back(hasTrackDetails ? cand->bestTrack()->dxy() : default_value);
            tauTuple().pfCand_dxy_error.push_back(hasTrackDetails ? cand->dxyError() : default_value);
            tauTuple().pfCand_dz.push_back(hasTrackDetails ? cand->bestTrack()->dz() : default_value);
            tauTuple().pfCand_dz_error.push_back(hasTrackDetails ? cand->dzError() : default_value);
            tauTuple().pfCand_track_chi2.push_back(
                        hasTrackDetails ? static_cast<float>(cand->bestTrack()->chi2()) : default_value);
            tauTuple().pfCand_track_ndof.push_back(
                        hasTrackDetails ? static_cast<float>(cand->bestTrack()->ndof()) : default_value);

            float hcal_fraction = cand->rawHcalEnergy()/(cand->rawHcalEnergy()+cand->rawEcalEnergy());
            tauTuple().pfCand_hcalFraction.push_back(hcal_fraction);
            tauTuple().pfCand_rawCaloFraction.push_back((cand->rawEcalEnergy()+cand->rawHcalEnergy())/cand->energy());

            //electron vars
            tauTuple().pfCand_ele_trackMomentumAtVtx.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->trackMomentumAtVtx().R() : default_value);
            tauTuple().pfCand_ele_trackMomentumAtCalo.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->trackMomentumAtCalo().R() : default_value);
            tauTuple().pfCand_ele_trackMomentumOut.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->trackMomentumOut().R() : default_value);
            tauTuple().pfCand_ele_trackMomentumAtEleClus.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->trackMomentumAtEleClus().R() : default_value);
            tauTuple().pfCand_ele_trackMomentumAtVtxWithConstraint.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->trackMomentumAtVtxWithConstraint().R() : default_value);
            tauTuple().pfCand_ele_ecalEnergy.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->ecalEnergy() : default_value);
            tauTuple().pfCand_ele_ecalEnergy_error.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->ecalEnergyError() : default_value);
            tauTuple().pfCand_ele_eSuperClusterOverP.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->eSuperClusterOverP() : default_value);
            tauTuple().pfCand_ele_eSeedClusterOverP.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->eSeedClusterOverP() : default_value);
            tauTuple().pfCand_ele_eSeedClusterOverPout.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->eSeedClusterOverPout() : default_value);
            tauTuple().pfCand_ele_eEleClusterOverPout.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->eEleClusterOverPout() : default_value);
            tauTuple().pfCand_ele_deltaEtaSuperClusterTrackAtVtx.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaEtaSuperClusterTrackAtVtx() : default_value);
            tauTuple().pfCand_ele_deltaEtaSeedClusterTrackAtCalo.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaEtaSeedClusterTrackAtCalo() : default_value);
            tauTuple().pfCand_ele_deltaEtaEleClusterTrackAtCalo.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaEtaEleClusterTrackAtCalo() : default_value);
            tauTuple().pfCand_ele_deltaPhiEleClusterTrackAtCalo.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaPhiEleClusterTrackAtCalo() : default_value);
            tauTuple().pfCand_ele_deltaPhiSuperClusterTrackAtVtx.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaPhiSuperClusterTrackAtVtx() : default_value);
            tauTuple().pfCand_ele_deltaPhiSeedClusterTrackAtCalo.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->deltaPhiSeedClusterTrackAtCalo() : default_value);

            tauTuple().pfCand_ele_mvaInput_earlyBrem.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->mvaInput().earlyBrem : default_value);
            tauTuple().pfCand_ele_mvaInput_lateBrem.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->mvaInput().lateBrem : default_value);
            tauTuple().pfCand_ele_mvaInput_sigmaEtaEta.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->mvaInput().sigmaEtaEta : default_value);
            tauTuple().pfCand_ele_mvaInput_hadEnergy.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->mvaInput().hadEnergy : default_value);
            tauTuple().pfCand_ele_mvaInput_deltaEta.push_back(cand->gsfElectronRef().isNonnull() ? cand->gsfElectronRef()->mvaInput().deltaEta : default_value);

            bool hasClosest = cand->gsfElectronRef().isNonnull() && cand->gsfElectronRef()->closestCtfTrackRef().isNonnull();
            tauTuple().pfCand_ele_closestCtfTrack_normalizedChi2.push_back(hasClosest
                        ? static_cast<float>(cand->gsfElectronRef()->closestCtfTrackRef()->normalizedChi2()) : default_value);
            tauTuple().pfCand_ele_closestCtfTrack_numberOfValidHits.push_back(
                        hasClosest ? cand->gsfElectronRef()->closestCtfTrackRef()->numberOfValidHits() : default_int_value);
        }
    }

    void FillElectrons(const std::vector<const reco::RecoEcalCandidate*>& electrons)
    {
        for(const reco::RecoEcalCandidate* ele : electrons) {
            tauTuple().ele_pt.push_back(static_cast<float>(ele->polarP4().pt()));
            tauTuple().ele_eta.push_back(static_cast<float>(ele->polarP4().eta()));
            tauTuple().ele_phi.push_back(static_cast<float>(ele->polarP4().phi()));
            tauTuple().ele_mass.push_back(static_cast<float>(ele->polarP4().mass()));
            float cc_ele_energy, cc_gamma_energy;
            int cc_n_gamma;
            CalculateElectronClusterVars(*ele, cc_ele_energy, cc_gamma_energy, cc_n_gamma);
            tauTuple().ele_cc_ele_energy.push_back(cc_ele_energy);
            tauTuple().ele_cc_gamma_energy.push_back(cc_gamma_energy);
            tauTuple().ele_cc_n_gamma.push_back(cc_n_gamma);

            const auto& gsfTrack = ele->gsfTrack();
            tauTuple().ele_gsfTrack_normalizedChi2.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->normalizedChi2()) : default_value);
            tauTuple().ele_gsfTrack_numberOfValidHits.push_back(
                        gsfTrack.isNonnull() ? gsfTrack->numberOfValidHits() : default_int_value);
            tauTuple().ele_gsfTrack_pt.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->pt()) : default_value);
            tauTuple().ele_gsfTrack_pt_error.push_back(
                        gsfTrack.isNonnull() ? static_cast<float>(gsfTrack->ptError()) : default_value);

        }
    }

    void FillMuons(const std::vector<const reco::Muon*>& muons, const std::vector<PFCandDesc>& cands)
    {
        for(const reco::Muon* muon : muons) {
            tauTuple().muon_pt.push_back(static_cast<float>(muon->polarP4().pt()));
            tauTuple().muon_eta.push_back(static_cast<float>(muon->polarP4().eta()));
            tauTuple().muon_phi.push_back(static_cast<float>(muon->polarP4().phi()));
            tauTuple().muon_mass.push_back(static_cast<float>(muon->polarP4().mass()));
            bool hasBestTrack = muon->bestTrack() != nullptr;
            tauTuple().muon_dxy.push_back(hasBestTrack ? static_cast<float>(muon->bestTrack()->dxy()) : default_value);
            tauTuple().muon_dxy_error.push_back(hasBestTrack ? static_cast<float>(std::abs(muon->bestTrack()->dxy())/muon->dxyError()) : default_value);

            const bool normalizedChi2_valid = hasBestTrack && muon->bestTrack()->outerOk() && muon->bestTrack()->normalizedChi2() >= 0;
            if(normalizedChi2_valid){
                tauTuple().muon_normalizedChi2.push_back(muon->bestTrack()->normalizedChi2());
                if(muon->bestTrack()->innerOk())
                    tauTuple().muon_numberOfValidHits.push_back(static_cast<int>(muon->bestTrack()->numberOfValidHits()));
            }

            double segmentCompatibility = muon::segmentCompatibility(*muon,reco::Muon::SegmentAndTrackArbitration);
            tauTuple().muon_segmentCompatibility.push_back(static_cast<float>(segmentCompatibility));
            tauTuple().muon_caloCompatibility.push_back(muon->caloCompatibility());

            bool pfEcalEnergy_valid = false;
            double pfEcalEnergy = 0.0;
            for (const PFCandDesc& cand_desc : cands) {
                const reco::PFCandidate* pfcand = cand_desc.candidate;
              if (pfcand->muonRef().isNonnull()) {
                if(muon == &(*pfcand->muonRef())){
                    pfEcalEnergy_valid = pfcand->ecalEnergy() >= 0;
                    pfEcalEnergy = pfcand->ecalEnergy();
                }
              }
            }

            tauTuple().muon_pfEcalEnergy.push_back(pfEcalEnergy_valid ? pfEcalEnergy : default_value);

            const MuonHitMatch hit_match(*muon);
            for(int subdet : MuonHitMatch::ConsideredSubdets()) {
                const std::string& subdetName = MuonHitMatch::SubdetName(subdet);
                for(int station = MuonHitMatch::first_station_id; station <= MuonHitMatch::last_station_id; ++station) {
                    const std::string matches_branch_name = "muon_n_matches_" + subdetName + "_"
                            + std::to_string(station);
                    const std::string hits_branch_name = "muon_n_hits_" + subdetName + "_" + std::to_string(station);

                    const unsigned n_matches = hit_match.NMatches(subdet, station);
                    const unsigned n_hits = hit_match.NHits(subdet, station);
                    tauTuple.get<std::vector<int>>(matches_branch_name).push_back(static_cast<int>(n_matches));
                    tauTuple.get<std::vector<int>>(hits_branch_name).push_back(static_cast<int>(n_hits));
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
    //edm::EDGetTokenT<pat::JetCollection> jets_token;
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

    TauTupleProducerData* data;
    tau_tuple::TauTuple& tauTuple;
    tau_tuple::SummaryTuple& summaryTuple;
};

} // namespace tau_analysis

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTupleProducer = tau_analysis::TauTupleProducer;
DEFINE_FWK_MODULE(TauTupleProducer);
