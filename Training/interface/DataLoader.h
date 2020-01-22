/*! Creates class to create vector for training.
*/

#include "TauML/Analysis/include/TrainTuple.h"
#include "TauML/Analysis/include/SummaryTuple.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/TextIO.h"


struct Data {
    Data(size_t n) : caloTowers(n), pixelTracks(n), isTau_vector(n) {}
    std::vector<float> caloTowers, pixelTracks;
    std::vector<bool> isTau;
};

struct Descriptor {
public:
    using Tau = train_tuple::Tau; //singolo evento
    using TrainTuple = train_tuple::TrainTuple;

    Descriptor(std::vector<std::string> _inputFile_names, bool _require_match) : {}

    const Event& getEvt() const { *tuple.get(); }
    const size_t getL1Index() const { currindex_l1_tau}
    bool GetNext()
    {
        while(GoToNext()) {
            if(event.l1Tau_pt.at(n) <= 20) continue;
            if(std::abs(event.l1Tau_eta.at(n)) >= 2.1) continue;
            if(require_match){
                for(unsigned h = 0; h < event.lepton_gen_match.size(); ++h){
                    if(event.lepton_gen_match.at(h) != 5) continue;
                    if(event.lepton_gen_visible_p4_pt.at(h) <= 20) continue;
                    if(std::abs(event.lepton_gen_visible_p4_eta.at(h)) >= 2.1) continue;
                    float delta_eta = event.l1Tau_eta.at(n) - event.lepton_gen_visible_p4_eta.at(h);
                    float delta_phi = ROOT::TVector2::Phi_mpi_pi(event.l1Tau_phi.at(n) - event.lepton_gen_visible_p4_phi.at(h));
                    float deltaR = std::sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
                    if(deltaR < 0.4){
                        matched_l1 = true;
                        return true;
                    }
                }
            }
            return true;
        }
        return false;

private:

    bool GoToNext()
    {
        if(!initialized){
            file_index = 0;
            index_tuple = 0;
            index_l1_tau = 0;
        }

        while(true) {
            if(index_l1_tau + 1 < (*ntuple)().l1Tau_pt.size()) {
                ++index_l1_tau;
                return true;
            }
            if(index_tuple < ntuple.GetEntries()) {
                ++index_tuple;
                ntuple.GetEntry(index_tuple);
                index_l1_tau = -1;
                continue;
            }
            if(file_index < inputFile_names.size()){
                ++file_index;
                inputFile_names.at(file_index);
                index_tuple = -1;
                index_l1_tau = -1;
                continue;
            }
            initialized = true;
            return false;
        }
    }

    std::vector<std::string> inputFile_names;
    bool require_match;
    size_t file_index;
    std::shared<TFile> file_opened;
    std::shared<TrainTuple> ntuple;
    Long64_t index_tuple;
    Int_t index_l1_tau;
    bool initialized;



};


class DataLoader {
public:
    using Tau = train_tuple::Tau; //singolo evento
    using TrainTuple = train_tuple::TrainTuple;

    DataLoader(std::vector<std::string> _inputFiles_signal, std::vector<std::string> _inputFiles_bkg, size_t _n) :
        n(_n)
    {
        file_id = rand()%2;
        if(file_id == 0)
            inputFiles = _inputFiles_signal;
        if(file_id == 1)
            inputFiles = _inputFiles_bkg;
    }

    Data GetData()
    {
        Data data(n);
        //Loop over 100K L1 taus
        for(unsigned f = 0; f < inputFiles.size(); ++f){
            auto originalFile = root_ext::OpenRootFile(inputFiles.at(f));
            TrainTuple trainTuple("taus", originalFile.get(), true);
            for(unsigned h = 0; h < n; ++h){
                size_t index_L1_tau = 0;
                Tau eventTau = GetTau(trainTuple,index_L1_tau);
                //loop on caloTower and fill the vector
                bool matched_calo = false;
                std::vector<float> calo_indexes;
                for(unsigned c = 0; c < eventTau.caloTower_energy.size(); ++c){
                    float delta_eta = eventTau.caloTower_eta.at(c) - eventTau.l1Tau_eta.at(index_L1_tau);
                    float delta_phi = ROOT::TVector2::Phi_mpi_pi(eventTau.l1Tau_phi.at(index_L1_tau) - eventTau.caloTower_phi.at(c));
                    float deltaR = std::sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
                    if(deltaR < 0.4){
                        calo_indexes.push_back(c);
                    }
                }

                const auto Comparitor_energy = [&](size_t h1, size_t h2) -> bool
                {
                    if(h1 == h2) return false;
                    float calo_energy_1 = eventTau.caloTower_energy.at(h1);
                    float calo_energy_2 = eventTau.caloTower_energy.at(h2);
                    return calo_energy_1 > calo_energy_2;
                };

                std::sort(calo_indexes.begin(),calo_indexes.end(),Comparitor_energy);
                calo_most_index = calo_indexes.at(0);
                for(unsigned index = 0; index < 10; ++index){
                    calo_index = calo_indexes.at(index);
                    data.caloTowers.push_back(eventTau.caloTower_eta.at(calo_most_index));
                    data.caloTowers.push_back(eventTau.caloTower_eta.at(calo_most_index)-eventTau.caloTower_eta.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_phi.at(calo_most_index)-eventTau.caloTower_phi.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_pt.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_eta.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_phi.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_energy.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_emEnergy.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_hadEnergy.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_outerEnergy.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_emPosition_x.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_emPosition_y.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_emPosition_z.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_hadPosition_x.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_hadPosition_y.at(calo_index));
                    data.caloTowers.push_back(eventTau.caloTower_hadPosition_z.at(calo_index));
                }
                //loop on pixelTracks and fill the vector
                std::vector<float> pixel_indexes;
                for(unsigned t = 0; t < eventTau.track_pt.size(); ++t){
                    float delta_eta = eventTau.track_eta.at(t) - eventTau.l1Tau_eta.at(index_L1_tau);
                    float delta_phi = ROOT::TVector2::Phi_mpi_pi(eventTau.l1Tau_phi.at(index_L1_tau) - eventTau.track_phi.at(t));
                    float deltaR = std::sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
                    if(deltaR < 0.4){
                        pixel_indexes.push_back(t);
                    }
                }

                const auto Comparitor_pt = [&](size_t h1, size_t h2) -> bool
                {
                    if(h1 == h2) return false;
                    float pt_1 = eventTau.track_pt.at(h1);
                    float pt_2 = eventTau.track_pt.at(h2);
                    return pt_1 > pt_2;
                };

                std::sort(pixel_indexes.begin(),pixel_indexes.end(),Comparitor_pt);
                for(unsigned id = 0; id < 10; ++id){
                    pixel_index = pixel_indexes.at(id);
                    data.pixelTracks.push_back(eventTau.caloTower_eta.at(calo_most_index));
                    data.pixelTracks.push_back(eventTau.caloTower_eta.at(calo_most_index)-eventTau.track_eta.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.caloTower_phi.at(calo_most_index)-eventTau.track_phi.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_pt.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_eta.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_phi.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_found.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_lost.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_chi2.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_ndof.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_charge.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_qualityMask.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_dxy.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_dz.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_vx.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_vy.at(pixel_index));
                    data.pixelTracks.push_back(eventTau.track_vz.at(pixel_index));
                }

                //isTau vector
                if(file_id == 0)
                    data.isTau.push_back(true);
                if(file_id == 1)
                    data.isTau.push_back(false);
            }

        }

        return data;
    }


    }



private:
    std::vector<std::string> inputFiles;
    size_t n; //max size of taus to be stored
    int file_id;
};
