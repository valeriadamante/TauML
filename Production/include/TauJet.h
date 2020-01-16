/* Tau jet candidate.
*/

#pragma once

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "GenTruthTools.h"

namespace tau_analysis {

enum class JetTauMatch { NoMatch = 0, PF = 1, dR = 2 };

struct PFCandDesc {
    const reco::PFCandidate* candidate;
    bool jetDaughter{false}, tauSignal{false}, leadChargedHadrCand{false}, tauIso{false};
};

struct TauJet {
    const pat::Jet* jet{nullptr};
    const reco::PFTau* tau{nullptr};
    JetTauMatch jetTauMatch{JetTauMatch::NoMatch};
    int jetIndex{-1}, tauIndex{-1};

    std::vector<PFCandDesc> cands;
    std::vector<const reco::RecoEcalCandidate*> electrons;
    std::vector<const reco::Muon*> muons;
    analysis::gen_truth::LeptonMatchResult jetGenLeptonMatchResult, tauGenLeptonMatchResult;
    analysis::gen_truth::QcdMatchResult jetGenQcdMatchResult, tauGenQcdMatchResult;

    TauJet(const pat::Jet* _jet, size_t _jetIndex);
    TauJet(const reco::PFTau* _tau, size_t _tauIndex);
    TauJet(const pat::Jet* _jet, const reco::PFTau* _tau, JetTauMatch _jetTauMatch, size_t _jetIndex, size_t _tauIndex);
};

struct TauJetBuilderSetup {
    double minJetPt{10.};
    double maxJetEta{3.};
    bool forceTauJetMatch{false}, useOnlyTauObjectMatch{false};
    double tauJetMatchDeltaR2Threshold{std::pow(0.2, 2)};
    double objectMatchDeltaR2ThresholdJet{std::pow(0.7, 2)};
    double objectMatchDeltaR2ThresholdTau{std::pow(0.5, 2)};
};

class TauJetBuilder {
public:
    using IndexSet = std::set<size_t>;
    using PolarLorentzVector = reco::LeafCandidate::PolarLorentzVector;

    TauJetBuilder(const TauJetBuilderSetup& setup, const pat::JetCollection& jets, const std::vector<reco::PFTau>& taus,
                  const std::vector<reco::PFCandidate>& cands, const std::vector<reco::RecoEcalCandidate>& electrons,
                  const reco::MuonCollection& muons, const reco::GenParticleCollection* genParticles);

    TauJetBuilder(const TauJetBuilder&) = delete;
    TauJetBuilder& operator=(const TauJetBuilder&) = delete;

    std::vector<TauJet> Build();

private:
    static bool IsJetDaughter(const pat::Jet& jet, const reco::PFCandidate& cand);
    static bool IsTauSignalCand(const reco::PFTau& tau, const reco::PFCandidate& cand);
    static bool IsTauIsoCand(const reco::PFTau& tau, const reco::PFCandidate& cand);
    static bool IsLeadChargedHadrCand(const reco::PFTau& tau, const reco::PFCandidate& cand);

    void MatchJetsAndTaus(JetTauMatch matchStrategy, std::vector<TauJet>& tauJets);
    bool FindJet(const reco::PFTau& tau, JetTauMatch matchStrategy, size_t& matchedJetIndex) const;
    bool GetMatchReferences(const pat::Jet* jet, const reco::PFTau* tau,
                            PolarLorentzVector& ref_p4, double& deltaR2) const;
    std::vector<PFCandDesc> FindMatchedPFCandidates(const pat::Jet* jet, const reco::PFTau* tau) const;
    std::vector<const reco::RecoEcalCandidate*> FindMatchedElectrons(const pat::Jet* jet, const reco::PFTau* tau) const;
    std::vector<const reco::Muon*> FindMatchedMuons(const pat::Jet* jet, const reco::PFTau* tau) const;

private:
    const TauJetBuilderSetup setup_;
    const pat::JetCollection& jets_;
    const std::vector<reco::PFTau>& taus_;
    const std::vector<reco::PFCandidate>& cands_;
    const std::vector<reco::RecoEcalCandidate>& electrons_;
    const reco::MuonCollection& muons_;
    const reco::GenParticleCollection* genParticles_;

    IndexSet availableJets_, processedJets_;
    IndexSet availableTaus_, processedTaus_;
};

} // namespace tau_analysis
