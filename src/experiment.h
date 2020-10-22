#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "mcts.h"
#include "simulator.h"
#include "statistic.h"
#include <fstream>
#include <string>
#include <vector>

////#ifdef _OPENMP
//#include<omp.h>
////#endif

//----------------------------------------------------------------------------

struct RESULTS
{
    void Clear();

    STATISTIC Time;
    STATISTIC Reward;
    STATISTIC DiscountedReturn;
    STATISTIC UndiscountedReturn;
    STATISTIC Cumulative_I;
    STATISTIC Total_Testing_Number;
    STATISTIC Reward_Numerator;
    STATISTIC Reward_Denomiantor;
    STATISTIC* S = NULL;
    STATISTIC* E = NULL;
    STATISTIC* I1 = NULL;
    STATISTIC* I2 = NULL;
    STATISTIC* H = NULL;
    STATISTIC* R = NULL;
    STATISTIC* testing_types = NULL;
    STATISTIC* testing_number = NULL;

    std::vector <double> all_S;
    std::vector <double> all_E;
    std::vector <double> all_I1;
    std::vector <double> all_I2;
    std::vector <double> all_H;
    std::vector <double> all_R;
    std::vector <double> all_testing_types;
    std::vector <double> all_testing_numbers;

    double initial_cum_S;
    double initial_cum_E;
    double initial_cum_I1;
    double initial_cum_I2;
    double initial_cum_H;
    double initial_cum_R;

};

inline void RESULTS::Clear()
{
    Time.Clear();
    Reward.Clear();
    DiscountedReturn.Clear();
    UndiscountedReturn.Clear();
    Cumulative_I.Clear();
    Total_Testing_Number.Clear();
    delete [] S;
    delete [] E;
    delete [] I1;
    delete [] I2;
    delete [] H;
    delete [] R;
    delete [] testing_types;
    delete [] testing_number;
}

//----------------------------------------------------------------------------

class EXPERIMENT
{
public:

    struct PARAMS
    {
        PARAMS();
        
        int NumRuns;
        int NumSteps;
        int SimSteps;
        double TimeOut;
        int MinDoubles, MaxDoubles;
        int TransformDoubles;
        int TransformAttempts;
        double Accuracy;
        int UndiscountedHorizon;
        bool AutoExploration;
        std::string problems_name;
        int fixed_action_ID;
        int is_random_action = false;
        int debug_mode;
    };

    EXPERIMENT(const SIMULATOR& real, const SIMULATOR& simulator, 
        const std::string& outputFile, 
        EXPERIMENT::PARAMS& expParams, MCTS::PARAMS& searchParams);

    void Run();
    void MultiRun();
    void DiscountedReturn();
    void AverageReward();

private:

    const SIMULATOR& Real;
    const SIMULATOR& Simulator;
    EXPERIMENT::PARAMS& ExpParams;
    MCTS::PARAMS& SearchParams;
    RESULTS Results;

    std::ofstream OutputFile;
    std::ofstream processOutputFile;
    std::ofstream allOutputFile;
};

//----------------------------------------------------------------------------

#endif // EXPERIMENT_H
