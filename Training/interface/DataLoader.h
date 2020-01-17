/*! Creates class to create vector for training.
*/

#include "TauML/Analysis/include/TrainTuple.h"
#include "TauML/Analysis/include/SummaryTuple.h"

class DataLoader {
public:
    using Tau = train_tuple::Tau;
    using TrainTuple = train_tuple::TrainTuple;

    DataLoader(std::vector<std::string> _inputFiles) : inputFiles(_inputFiles) { }

    void GetData()
    {
        //Loop over 100K L1 taus
    }

    size_t GetTau()
    {
        //Loop over events
        //return index of L1 tau and event
    }



private:
    std::vector<std::string> inputFiles;
};
