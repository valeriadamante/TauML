/*! Definition of a tuple with all event information that is required for the tau analysis.
*/

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "TauML/Analysis/include/TauIdResults.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include <Math/VectorUtil.h>

#define TAU_ID(name, pattern, has_raw, wp_list) VAR(std::vector<uint16_t>, name) VAR(std::vector<Float_t>, name##raw)
#define TAU_VAR(type, name) VAR(std::vector<type>, tau_##name)
#define CAND_VAR(type, name) VAR(std::vector<type>, pfCand_##name)
#define ELE_VAR(type, name) VAR(std::vector<type>, ele_##name)
#define MUON_VAR(type, name) VAR(std::vector<type>, muon_##name)
#define CALO_TOWER_VAR(type, name) VAR(std::vector<type>, caloTower_##name)
#define CALO_TAU_VAR(type, name) VAR(std::vector<type>, caloTau_##name)
#define TRACK_VAR(type, name) VAR(std::vector<type>, track_##name)

#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define CAND_VAR2(type, name1, name2) CAND_VAR(type, name1) CAND_VAR(type, name2)
#define CAND_VAR3(type, name1, name2, name3) CAND_VAR2(type, name1, name2) CAND_VAR(type, name3)
#define CAND_VAR4(type, name1, name2, name3, name4) CAND_VAR3(type, name1, name2, name3) CAND_VAR(type, name4)

#define ELE_VAR2(type, name1, name2) ELE_VAR(type, name1) ELE_VAR(type, name2)
#define ELE_VAR3(type, name1, name2, name3) ELE_VAR2(type, name1, name2) ELE_VAR(type, name3)
#define ELE_VAR4(type, name1, name2, name3, name4) ELE_VAR3(type, name1, name2, name3) ELE_VAR(type, name4)

#define MUON_VAR2(type, name1, name2) MUON_VAR(type, name1) MUON_VAR(type, name2)
#define MUON_VAR3(type, name1, name2, name3) MUON_VAR2(type, name1, name2) MUON_VAR(type, name3)
#define MUON_VAR4(type, name1, name2, name3, name4) MUON_VAR3(type, name1, name2, name3) MUON_VAR(type, name4)

#define TAU_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, npv) /* number of primary vertices */ \
    VAR(Float_t, rho) /* fixed grid energy density */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, trainingWeight) /* training weight */ \
    VAR(Int_t, sampleType) /* type of the sample (MC, Embedded or Data) */ \
    VAR(Float_t, npu) /* number of in-time pu interactions added to the event */ \
    VAR3(Float_t, pv_x, pv_y, pv_z) /* position of the primary vertex (PV) */ \
    VAR(Float_t, pv_chi2) /* chi^2 of the primary vertex (PV) */ \
    VAR(Float_t, pv_ndof) /* number of degrees of freedom of the primary vertex (PV) */ \
    /* Basic tau variables */ \
    TAU_VAR(Int_t, index) /* index of the tau */ \
    TAU_VAR(Float_t, pt) /* 4-momentum of the tau */ \
    TAU_VAR(Float_t, eta) /* 4-momentum of the tau */ \
    TAU_VAR(Float_t, phi) /* 4-momentum of the tau */ \
    TAU_VAR(Float_t, mass) /* 4-momentum of the tau */ \
    TAU_VAR(Int_t, charge) /* tau charge */ \
    /* Lepton gen match variables */ \
    VAR(std::vector<Int_t>, lepton_gen_match) /* matching with leptons on the generator level */\
    VAR(std::vector<Int_t>, lepton_genLast_charge) /* charge of the matched gen last copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genLast_pt) /* 4-momentum of the matched gen last copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genLast_eta) /* 4-momentum of the matched gen last copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genLast_phi) /* 4-momentum of the matched gen last copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genLast_mass) /* 4-momentum of the matched gen last copy lepton */ \
    VAR(std::vector<Int_t>, lepton_genFirst_charge) /* charge of the matched gen first copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genFirst_pt) /* 4-momentum of the matched gen first copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genFirst_eta) /* 4-momentum of the matched gen first copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genFirst_phi) /* 4-momentum of the matched gen first copy lepton */ \
    VAR(std::vector<Float_t>, lepton_genFirst_mass) /* 4-momentum of the matched gen first copy lepton */ \
    VAR(std::vector<Int_t>, lepton_gen_vis_daug_index) /* index of the matched lepton mother */ \
    VAR(std::vector<Int_t>, lepton_gen_vis_daug_pdg) /* PDG of the matched lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_daug_pt) /* 4-momenta of the visible products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_daug_eta) /* 4-momenta of the visible products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_daug_phi) /* 4-momenta of the visible products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_daug_mass) /* 4-momenta of the visible products of the matched gen lepton */ \
    VAR(std::vector<Int_t>, lepton_gen_vis_rad_index) /* index of the matched lepton mother */ \
    VAR(std::vector<Int_t>, lepton_gen_vis_rad_pdg) /* PDG of the matched lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_rad_pt) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_rad_eta) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_rad_phi) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_vis_rad_mass) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_p4_pt) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_p4_eta) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_p4_phi) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_p4_mass) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_rad_p4_pt) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_rad_p4_eta) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_rad_p4_phi) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<Float_t>, lepton_gen_visible_rad_p4_mass) /* 4-momenta of the visible rad products of the matched gen lepton */ \
    VAR(std::vector<UInt_t>, lepton_gen_n_charged_hadrons) /* n_charged_hadrons of the matched lepton */ \
    VAR(std::vector<UInt_t>, lepton_gen_n_neutral_hadrons) /* n_neutral_hadrons of the matched lepton */ \
    VAR(std::vector<UInt_t>, lepton_gen_n_gammas) /* n_gammas of the matched lepton */ \
    VAR(std::vector<UInt_t>, lepton_gen_n_gammas_rad) /* n_gammas_rad of the matched lepton */ \
    /* Tau ID variables */ \
    TAU_VAR(Int_t, decayMode) /* tau decay mode */ \
    TAU_VAR(Int_t, decayModeFinding) /* tau passed the old decay mode finding requirements */ \
    TAU_VAR(Int_t, decayModeFindingNewDMs) /* tau passed the new decay mode finding requirements */ \
    TAU_VAR(Float_t, chargedIsoPtSum) /* sum of the transverse momentums of charged pf candidates inside the tau isolation cone with dR < 0.5 */ \
    TAU_VAR(Float_t, chargedIsoPtSumdR03) /* sum of the transverse momentums of charged pf candidates inside the tau isolation cone with dR < 0.3 */ \
    TAU_VAR(Float_t, footprintCorrection) /* tau footprint correction inside the tau isolation cone with dR < 0.5 */ \
    TAU_VAR(Float_t, footprintCorrectiondR03) /* tau footprint correction inside the tau isolation cone with dR < 0.3 */ \
    TAU_VAR(Float_t, neutralIsoPtSum) /* sum of the transverse momentums of neutral pf candidates inside
                                     the tau isolation cone with dR < 0.5 */ \
    TAU_VAR(Float_t, neutralIsoPtSumWeight) /* weighted sum of the transverse momentums of neutral pf candidates inside
                                           the tau isolation cone with dR < 0.5 */ \
    TAU_VAR(Float_t, neutralIsoPtSumWeightdR03) /* weighted sum of the transverse momentums of neutral pf candidates inside
                                               the tau isolation cone with dR < 0.3 */ \
    TAU_VAR(Float_t, neutralIsoPtSumdR03) /* sum of the transverse momentums of neutral pf candidates inside
                                         the tau isolation cone with dR < 0.3 */ \
    TAU_VAR(Float_t, photonPtSumOutsideSignalCone) /* sum of the transverse momentums of photons
                                                  inside the tau isolation cone with dR < 0.5 */ \
    TAU_VAR(Float_t, photonPtSumOutsideSignalConedR03) /* sum of the transverse momentums of photons inside
                                                      the tau isolation cone with dR < 0.3 */ \
    TAU_VAR(Float_t, puCorrPtSum) /* pile-up correction for the sum of the transverse momentums */ \
    TAU_IDS() \
    /* Tau transverse impact paramters.
       See cmssw/RecoTauTag/RecoTau/plugins/PFTauTransverseImpactParameters.cc for details */ \
    /* VAR3(Float_t, tau_dxy_pca_x, tau_dxy_pca_y, tau_dxy_pca_z) The point of closest approach (PCA) of the leadPFChargedHadrCand to the primary vertex */ \
    TAU_VAR(Float_t, dxy) /* tau signed transverse impact parameter wrt to the primary vertex */ \
    TAU_VAR(Float_t, dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    TAU_VAR(Float_t, ip3d) /* tau signed 3D impact parameter wrt to the primary vertex */ \
    TAU_VAR(Float_t, ip3d_error) /* uncertainty of the 3D impact parameter measurement */ \
    TAU_VAR(Float_t, dz) /* tau dz of the leadChargedHadrCand wrt to the primary vertex */ \
    TAU_VAR(Float_t, dz_error) /* uncertainty of the tau dz measurement */ \
    TAU_VAR(Int_t, hasSecondaryVertex) /* tau has the secondary vertex */ \
    TAU_VAR(Float_t, sv_x) /* position of the secondary vertex */ \
    TAU_VAR(Float_t, sv_y) /* position of the secondary vertex */ \
    TAU_VAR(Float_t, sv_z) /* position of the secondary vertex */ \
    TAU_VAR(Float_t, flightLength_x) /* flight length of the tau */ \
    TAU_VAR(Float_t, flightLength_y) /* flight length of the tau */ \
    TAU_VAR(Float_t, flightLength_z) /* flight length of the tau */ \
    TAU_VAR(Float_t, flightLength_sig) /* significance of the flight length measurement */ \
    /* Extended tau variables */ \
    TAU_VAR(Float_t, pt_weighted_deta_strip) /* sum of pt weighted values of deta relative to tau candidate
                                                for all pf photon candidates, which are associated to signal */ \
    TAU_VAR(Float_t, pt_weighted_dphi_strip) /* sum of pt weighted values of dphi relative to tau candidate
                                                for all pf photon candidates, which are associated to signal */ \
    TAU_VAR(Float_t, pt_weighted_dr_signal) /* sum of pt weighted values of dr relative to tau candidate
                                               for all pf photon candidates, which are associated to signal */ \
    TAU_VAR(Float_t, pt_weighted_dr_iso) /* sum of pt weighted values of dr relative to tau candidate
                                            for all pf photon candidates, which are inside an isolation cone
                                            but not associated to signal */ \
    TAU_VAR(Float_t, leadingTrackNormChi2) /* normalized chi2 of leading track */ \
    TAU_VAR(Float_t, e_ratio) /* ratio of energy in ECAL over sum of energy in ECAL and HCAL */ \
    TAU_VAR(Float_t, gj_angle_diff) /* Gottfried-Jackson angle difference
                                       (defined olny when the secondary vertex is reconstructed) */ \
    TAU_VAR(Int_t, n_photons) /* total number of pf photon candidates with pT>500 MeV,
                                 which are associated to signal */ \
    TAU_VAR(Float_t, emFraction) /* tau->emFraction_MVA */ \
    TAU_VAR(Int_t, inside_ecal_crack) /* tau is inside the ECAL crack (1.46 < |eta| < 1.558) */ \
    TAU_VAR(Float_t, leadChargedCand_etaAtEcalEntrance) /* eta at ECAL entrance of the leadChargedCand */ \
    /* L1 objects */ \
    VAR(std::vector<Float_t>, l1Tau_pt) /* L1 pt candidate*/ \
    VAR(std::vector<Float_t>, l1Tau_eta) /* L1 eta candidate*/ \
    VAR(std::vector<Float_t>, l1Tau_phi) /* L1 phi candidate*/ \
    VAR(std::vector<Float_t>, l1Tau_mass) /* L1 mass candidate*/ \
    VAR(std::vector<Float_t>, l1Tau_hwIso) /* L1 hwIso candidate*/ \
    VAR(std::vector<Float_t>, l1Tau_hwQual) /* L1 quality candidate*/ \
    /* caloTower candidates */ \
    CALO_TOWER_VAR(Float_t, pt) /* caloTower pt candidate*/ \
    CALO_TOWER_VAR(Float_t, eta) /* caloTower eta candidate*/ \
    CALO_TOWER_VAR(Float_t, phi) /* caloTower phi candidate*/ \
    CALO_TOWER_VAR(Float_t, energy) /* caloTower energy candidate*/ \
    CALO_TOWER_VAR(Float_t, emEnergy) /* caloTower emEnergy candidate*/ \
    CALO_TOWER_VAR(Float_t, hadEnergy) /* caloTower hadEnergy candidate*/ \
    CALO_TOWER_VAR(Float_t, outerEnergy) /* caloTower outerEnergy candidate*/ \
    CALO_TOWER_VAR(Float_t, emPosition_x) /* caloTower emPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, emPosition_y) /* caloTower emPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, emPosition_z) /* caloTower emPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, hadPosition_x) /* caloTower hadPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, hadPosition_y) /* caloTower hadPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, hadPosition_z) /* caloTower hadPosition candidate*/ \
    CALO_TOWER_VAR(Float_t, hadEnergyHeOuterLayer) /* caloTower hadEnergyHeOuterLayer candidate*/ \
    CALO_TOWER_VAR(Float_t, hadEnergyHeInnerLayer) /* caloTower hadEnergyHeInnerLayer candidate*/ \
    CALO_TOWER_VAR(Float_t, energyInHB) /* caloTower energyInHB candidate*/ \
    CALO_TOWER_VAR(Float_t, energyInHE) /* caloTower energyInHE candidate*/ \
    CALO_TOWER_VAR(Float_t, energyInHF) /* caloTower energyInHF candidate*/ \
    CALO_TOWER_VAR(Float_t, energyInHO) /* caloTower energyInHO candidate*/ \
    CALO_TOWER_VAR(Int_t, numBadEcalCells) /* caloTower numBadEcalCells candidate*/ \
    CALO_TOWER_VAR(Int_t, numRecoveredEcalCells) /* caloTower numRecoveredEcalCells candidate*/ \
    CALO_TOWER_VAR(Int_t, numProblematicEcalCells) /* caloTower numProblematicEcalCells candidate*/ \
    CALO_TOWER_VAR(Int_t, numBadHcalCells) /* caloTower numBadHcalCells candidate*/ \
    CALO_TOWER_VAR(Int_t, numRecoveredHcalCells) /* caloTower numRecoveredHcalCells candidate*/ \
    CALO_TOWER_VAR(Int_t, numProblematicHcalCells) /* caloTower numProblematicHcalCells candidate*/ \
    CALO_TOWER_VAR(Float_t, ecalTime) /* caloTower ecalTime candidate*/ \
    CALO_TOWER_VAR(Float_t, hcalTime) /* caloTower hcalTime candidate*/ \
    /* pixelTracks candidates */ \
    TRACK_VAR(Float_t, pt) /* pixeltrack pt candidate*/ \
    TRACK_VAR(Float_t, eta) /* pixeltrack eta candidate*/ \
    TRACK_VAR(Float_t, phi) /* pixeltrack phi candidate*/ \
    TRACK_VAR(Bool_t, outerOk) /* pixeltrack outerOk candidate*/ \
    TRACK_VAR(Bool_t, innerOk) /* pixeltrack innerOk candidate*/ \
    TRACK_VAR(Int_t, found) /* pixeltrack Number of valid hits on track candidate*/ \
    TRACK_VAR(Int_t, lost) /* pixeltrack Number of lost (=invalid) hits on track candidate*/ \
    TRACK_VAR(Float_t, chi2) /* pixeltrack chi2 candidate*/ \
    TRACK_VAR(Float_t, ndof) /* pixeltrack ndof candidate*/ \
    TRACK_VAR(Int_t, charge) /* pixelTrack charge candidate*/ \
    TRACK_VAR(Int_t, algo) /* pixelTrack algo candidate*/ \
    TRACK_VAR(uint16_t, qualityMask) /* pixelTrack qualityMask candidate*/ \
    TRACK_VAR(Float_t, dxy) /* pixeltrack dxy candidate*/ \
    TRACK_VAR(Float_t, dz) /* pixeltrack dz candidate*/ \
    TRACK_VAR(Float_t, vx) /* pixeltrack vx candidate*/ \
    TRACK_VAR(Float_t, vy) /* pixeltrack vy candidate*/ \
    TRACK_VAR(Float_t, vz) /* pixeltrack vz candidate*/ \
    TRACK_VAR(Float_t, ptError) /* pixeltrack ptError candidate*/ \
    TRACK_VAR(Float_t, etaError) /* pixeltrack etaError candidate*/ \
    TRACK_VAR(Float_t, phiError) /* pixeltrack phiError candidate*/ \
    TRACK_VAR(Float_t, dxyError) /* pixeltrack dxyError candidate*/ \
    TRACK_VAR(Float_t, dzError) /* pixeltrack dzError candidate*/ \
    /* CaloTaus candidates */ \
    CALO_TAU_VAR(Float_t, pt) /* caloTau pt candidate*/ \
    CALO_TAU_VAR(Float_t, eta) /* caloTau eta candidate*/ \
    CALO_TAU_VAR(Float_t, phi) /* caloTau phi candidate*/ \
    CALO_TAU_VAR(Float_t, energy) /* caloTau energy candidate*/ \
    CALO_TAU_VAR(Float_t, maxEInEmTowers) /* caloTau maximum energy deposited in ECAL towers*/ \
    CALO_TAU_VAR(Float_t, maxEInHadTowers) /* caloTau maximum energy deposited in HCAL towers*/ \
    CALO_TAU_VAR(Float_t, energyFractionHadronic) /* caloTau jet hadronic energy fraction*/ \
    CALO_TAU_VAR(Float_t, emEnergyFraction) /* caloTau jet electromagnetic energy fraction*/ \
    CALO_TAU_VAR(Float_t, hadEnergyInHB) /* caloTau jet hadronic energy in HB*/ \
    CALO_TAU_VAR(Float_t, hadEnergyInHO) /* caloTau jet hadronic energy in HO*/ \
    CALO_TAU_VAR(Float_t, hadEnergyInHE) /* caloTau jet hadronic energy in HE*/ \
    CALO_TAU_VAR(Float_t, hadEnergyInHF) /* caloTau jet hadronic energy in HF*/ \
    CALO_TAU_VAR(Float_t, emEnergyInEB) /* caloTau jet electromagnetic energy in EB*/ \
    CALO_TAU_VAR(Float_t, emEnergyInEE) /* caloTau jet electromagnetic energy in EE*/ \
    CALO_TAU_VAR(Float_t, emEnergyInHF) /* caloTau jet electromagnetic energy extracted from HF*/ \
    CALO_TAU_VAR(Float_t, towersArea) /* caloTau area of contributing towers*/ \
    CALO_TAU_VAR(Int_t, n90) /* caloTau number of constituents carrying a 90% of the total Jet energy*/ \
    CALO_TAU_VAR(Int_t, n60) /* caloTau number of constituents carrying a 60% of the total Jet energy*/ \
    /* PF candidates */ \
    CAND_VAR4(Float_t, pt, eta, phi, mass) /* 4-momentum of the PF candidate */ \
    CAND_VAR(Int_t, pvAssociationQuality) /* information about how the association to the PV is obtained:
                                             NotReconstructedPrimary = 0, OtherDeltaZ = 1, CompatibilityBTag = 4,
                                             CompatibilityDz = 5, UsedInFitLoose = 6, UsedInFitTight = 7 */ \
    CAND_VAR(Int_t, fromPV) /* the association to PV=ipv. >=PVLoose corresponds to JME definition,
                               >=PVTight to isolation definition:
                               NoPV = 0, PVLoose = 1, PVTight = 2, PVUsedInFit = 3 */ \
    CAND_VAR(Float_t, puppiWeight) /* weight from full PUPPI */ \
    CAND_VAR(Float_t, puppiWeightNoLep) /* weight from PUPPI removing leptons */ \
    CAND_VAR(Int_t, pdgId) /* PDG identifier */ \
    CAND_VAR(Int_t, charge) /* electric charge */ \
    CAND_VAR(Int_t, lostInnerHits) /* enumerator specifying the number of lost inner hits:
                                      validHitInFirstPixelBarrelLayer = -1, noLostInnerHits = 0 (it could still not
                                      have a hit in the first layer, e.g. if it crosses an inactive sensor),
                                      oneLostInnerHit = 1, moreLostInnerHits = 2 */ \
    CAND_VAR(Int_t, numberOfPixelHits) /* number of valid pixel hits */ \
    CAND_VAR3(Float_t, vertex_x, vertex_y, vertex_z) /* position of the vertex to which the candidate is associated */ \
    CAND_VAR(Int_t, hasTrackDetails) /* has track details */ \
    CAND_VAR(Float_t, dxy) /* signed transverse impact parameter wrt to the primary vertex */ \
    CAND_VAR(Float_t, dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    CAND_VAR(Float_t, dz) /* dz wrt to the primary vertex */ \
    CAND_VAR(Float_t, dz_error) /* uncertainty of the dz measurement */ \
    CAND_VAR(Float_t, track_chi2) /* chi^2 of the pseudo track made with the candidate kinematics */ \
    CAND_VAR(Float_t, track_ndof) /* number of degrees of freedom of the pseudo track
                                     made with the candidate kinematics */ \
    CAND_VAR(Float_t, hcalFraction) /* fraction of ECAL and HCAL for HF and neutral hadrons
                                       and isolated charged hadrons */ \
    CAND_VAR(Float_t, rawCaloFraction) /* raw ECAL+HCAL energy over candidate energy for isolated charged hadrons */ \
    /* PAT electrons + pfCand ele */ \
    ELE_VAR4(Float_t, pt, eta, phi, mass) /* 4-momentum of the electron */ \
    ELE_VAR(Float_t, cc_ele_energy) /* energy of the first calo cluster in the electron super cluster */ \
    ELE_VAR(Float_t, cc_gamma_energy) /* sum of the energies of additional calo clusters
                                         in the electron super cluster */ \
    ELE_VAR(Int_t, cc_n_gamma) /* number of additional calo clusters in the electron super cluster */ \
    CAND_VAR(Float_t, ele_trackMomentumAtVtx) /* module of the track momentum at the PCA to the beam spot */ \
    CAND_VAR(Float_t, ele_trackMomentumAtCalo) /* module of the track momentum extrapolated at the supercluster position
                                             from the innermost track state */ \
    CAND_VAR(Float_t, ele_trackMomentumOut) /* module of the track momentum extrapolated at the seed cluster position
                                          from the outermost track state */ \
    CAND_VAR(Float_t, ele_trackMomentumAtEleClus) /* module of the track momentum extrapolated at the ele cluster position
                                                from the outermost track state */ \
    CAND_VAR(Float_t, ele_trackMomentumAtVtxWithConstraint) /* module of the track momentum at the PCA to the beam spot
                                                          using bs constraint */ \
    CAND_VAR(Float_t, ele_ecalEnergy) /*  corrected ECAL energy */ \
    CAND_VAR(Float_t, ele_ecalEnergy_error) /* uncertanty of the ECAL energy measurement */ \
    CAND_VAR(Float_t, ele_eSuperClusterOverP) /* supercluster energy / track momentum at the PCA to the beam spot */ \
    CAND_VAR(Float_t, ele_eSeedClusterOverP) /* seed cluster energy / track momentum at the PCA to the beam spot */ \
    CAND_VAR(Float_t, ele_eSeedClusterOverPout) /* seed cluster energy / track momentum at calo extrapolated
                                              from the outermost track state */ \
    CAND_VAR(Float_t, ele_eEleClusterOverPout) /* electron cluster energy / track momentum at calo extrapolated
                                             from the outermost track state */ \
    CAND_VAR(Float_t, ele_deltaEtaSuperClusterTrackAtVtx) /* supercluster eta - track eta position at calo extrapolated
                                                        from innermost track state */ \
    CAND_VAR(Float_t, ele_deltaEtaSeedClusterTrackAtCalo) /* seed cluster eta - track eta position at calo extrapolated
                                                        from the outermost track state */ \
    CAND_VAR(Float_t, ele_deltaEtaEleClusterTrackAtCalo) /* electron cluster eta - track eta position at calo extrapolated
                                                       from the outermost state */ \
    CAND_VAR(Float_t, ele_deltaPhiEleClusterTrackAtCalo) /* electron cluster phi - track phi position at calo extrapolated
                                                       from the outermost track state */ \
    CAND_VAR(Float_t, ele_deltaPhiSuperClusterTrackAtVtx) /* supercluster phi - track phi position at calo extrapolated
                                                        from the innermost track state */ \
    CAND_VAR(Float_t, ele_deltaPhiSeedClusterTrackAtCalo) /* seed cluster phi - track phi position at calo extrapolated
                                                        from the outermost track state */ \
    CAND_VAR(Int_t, ele_mvaInput_earlyBrem) /* early bremsstrahlung is detected:
                                                              unknown = -2, could not be evaluated = -1,
                                                              wrong = 0, true = 1 */ \
    CAND_VAR(Int_t, ele_mvaInput_lateBrem) /* late bremsstrahlung is detected:
                                                            unknown = -2, could not be evaluated = -1,
                                                            wrong = 0, true = 1 */ \
    CAND_VAR(Float_t, ele_mvaInput_sigmaEtaEta) /* Sigma-eta-eta with the PF cluster */ \
    CAND_VAR(Float_t, ele_mvaInput_hadEnergy) /* Associated PF Had Cluster energy */ \
    CAND_VAR(Float_t, ele_mvaInput_deltaEta) /* PF-cluster GSF track delta-eta */ \
    ELE_VAR(Float_t, gsfTrack_normalizedChi2) /* chi^2 divided by number of degrees of freedom of the GSF track */ \
    ELE_VAR(Int_t, gsfTrack_numberOfValidHits) /* number of valid hits on the GSF track */ \
    ELE_VAR(Float_t, gsfTrack_pt) /* pt of the GSF track */ \
    ELE_VAR(Float_t, gsfTrack_pt_error) /* uncertainty of the pt measurement of the GSF track */ \
    CAND_VAR(Float_t, ele_closestCtfTrack_normalizedChi2) /* chi^2 divided by number of degrees of freedom
                                                        of the closest CTF track */ \
    CAND_VAR(Int_t, ele_closestCtfTrack_numberOfValidHits) /* number of valid hits on the closest CTF track */ \
    /* PAT muons */ \
    MUON_VAR4(Float_t, pt, eta, phi, mass) /* 4-momentum of the muon */ \
    MUON_VAR(Float_t, dxy) /* signed transverse impact parameter of the inner track wrt to the primary vertex */ \
    MUON_VAR(Float_t, dxy_error) /* uncertainty of the transverse impact parameter measurement */ \
    MUON_VAR(Float_t, normalizedChi2) /* chi^2 divided by number of degrees of freedom of the global track */ \
    MUON_VAR(Int_t, numberOfValidHits) /* number of valid hits on the global track */ \
    MUON_VAR(Float_t, segmentCompatibility) /* segment compatibility for a track with matched muon info */ \
    MUON_VAR(Float_t, caloCompatibility) /* relative likelihood based on ECAL, HCAL, HO energy defined as
                                            L_muon / (L_muon + L_not_muon) */ \
    MUON_VAR(Float_t, pfEcalEnergy) /* PF based energy deposition in the ECAL */ \
    MUON_VAR4(Int_t, n_matches_DT_1, n_matches_DT_2, n_matches_DT_3, \
                     n_matches_DT_4) /* number of segment matches for the DT subdetector stations */ \
    MUON_VAR4(Int_t, n_matches_CSC_1, n_matches_CSC_2, n_matches_CSC_3, \
                     n_matches_CSC_4) /* number of segment matches for the CSC subdetector stations */ \
    MUON_VAR4(Int_t, n_matches_RPC_1, n_matches_RPC_2, n_matches_RPC_3, \
                     n_matches_RPC_4) /* number of segment matches for the RPC subdetector stations */ \
    MUON_VAR4(Int_t, n_hits_DT_1, n_hits_DT_2, n_hits_DT_3, \
                     n_hits_DT_4) /* number of valid and bad hits for the DT subdetector stations */ \
    MUON_VAR4(Int_t, n_hits_CSC_1, n_hits_CSC_2, n_hits_CSC_3, \
                     n_hits_CSC_4) /* number of valid and bad hits for the CSC subdetector stations */ \
    MUON_VAR4(Int_t, n_hits_RPC_1, n_hits_RPC_2, n_hits_RPC_3, \
                     n_hits_RPC_4) /* number of valid and bad hits for the RPC subdetector stations */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(train_tuple, Tau, TrainTuple, TAU_DATA, "taus")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(train_tuple, TrainTuple, TAU_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef TAU_DATA
#undef CAND_VAR
#undef CAND_VAR2
#undef CAND_VAR3
#undef CAND_VAR4
#undef ELE_VAR
#undef ELE_VAR2
#undef ELE_VAR3
#undef ELE_VAR4
#undef MUON_VAR
#undef MUON_VAR2
#undef MUON_VAR3
#undef MUON_VAR4
#undef TAU_ID
#undef TAU_VAR
#undef CALO_TOWER_VAR
#undef CALO_TAU_VAR
#undef TRACK_VAR

namespace train_tuple {

template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }
template<>
constexpr float DefaultFillValue<float>() { return -999.; }
template<>
constexpr int DefaultFillValue<int>() { return -999; }
template<>
constexpr unsigned DefaultFillValue<unsigned>() { return 0; }

} // namespace train_tuple
