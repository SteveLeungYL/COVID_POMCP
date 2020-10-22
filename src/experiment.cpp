#include "experiment.h"
#include "boost/timer.hpp"
#include "pomcp.h"
#include <random>

using namespace std;

EXPERIMENT::PARAMS::PARAMS()
        :   NumRuns(1000),
            NumSteps(100000),
            SimSteps(1000),
            TimeOut(3600),
            MinDoubles(0),
            MaxDoubles(20),
            TransformDoubles(-4),
            TransformAttempts(1000),
            Accuracy(0.01),
            UndiscountedHorizon(1000),
            AutoExploration(true)
{
}

EXPERIMENT::EXPERIMENT(const SIMULATOR& real,
                       const SIMULATOR& simulator, const string& outputFile,
                       EXPERIMENT::PARAMS& expParams, MCTS::PARAMS& searchParams)
        :   Real(real),
            Simulator(simulator),
            OutputFile(outputFile.c_str()),
            ExpParams(expParams),
            SearchParams(searchParams)
{
//// User added code:
    if (expParams.problems_name == "covid") {
        processOutputFile.open(outputFile + ".process.csv");
        allOutputFile.open(outputFile + ".process_all.csv");
    }

    if (ExpParams.AutoExploration)
    {
        if (SearchParams.UseRave)
            SearchParams.ExplorationConstant = 0;
        else
            SearchParams.ExplorationConstant = simulator.GetRewardRange();
    }
    MCTS::InitFastUCB(SearchParams.ExplorationConstant);
}

void EXPERIMENT::Run()
{
    boost::timer timer;

    MCTS mcts(Simulator, SearchParams);

    double undiscountedReturn = 0.0;
    double discountedReturn = 0.0;
    double discount = 1.0;
    bool termination = false;
    bool outOfParticles = false;
    int t = 0;

    int total_testing_number = 0;

    STATE* state = Real.CreateStartState();
    if (SearchParams.Verbose >= 1)
        Real.DisplayState(*state, cout);

    POMCP_STATE& pomcp_state = safe_cast<POMCP_STATE&>(*state);
    this->Results.initial_cum_S = pomcp_state.S;
    this->Results.initial_cum_E = pomcp_state.E;
    this->Results.initial_cum_I1 = pomcp_state.I1;
    this->Results.initial_cum_I2 = pomcp_state.I2;
    this->Results.initial_cum_H = pomcp_state.H;
    this->Results.initial_cum_R = pomcp_state.R;

    for (t = 0; t < ExpParams.NumSteps; t++)
    {
        int observation;
        double reward;
        int action = mcts.SelectAction();

        if (ExpParams.fixed_action_ID != 0) action = ExpParams.fixed_action_ID; // fix_action_ID = 0, means not limited.
        if (ExpParams.is_random_action){
            random_device rd;
            mt19937_64 gen(rd());
            uniform_real_distribution<float> rand_dis(0.0, 99.1);
            action = (int)rand_dis(gen);
        }

        termination = Real.Step(*state, action, observation, reward, cout);

        Results.Reward.Add(reward);
        undiscountedReturn += reward;
        discountedReturn += reward * discount;
        discount *= Real.GetDiscount();

        // Self implemented code:
        if (this->ExpParams.problems_name == "covid") {

            POMCP_STATE& pomcp_state = safe_cast<POMCP_STATE&>(*state);
            this->Results.S[t].Add(pomcp_state.S);
            this->Results.E[t].Add(pomcp_state.E);
            this->Results.I1[t].Add(pomcp_state.I1);
            this->Results.I2[t].Add(pomcp_state.I2);
            this->Results.H[t].Add(pomcp_state.H);
            this->Results.R[t].Add(pomcp_state.R);

            this->Results.all_S.push_back(pomcp_state.S);
            this->Results.all_E.push_back(pomcp_state.E);
            this->Results.all_I1.push_back(pomcp_state.I1);
            this->Results.all_I2.push_back(pomcp_state.I2);
            this->Results.all_H.push_back(pomcp_state.H);
            this->Results.all_R.push_back(pomcp_state.R);

            double testing_types = 0;
            int testing_number = 0;

            Real.Get_testing_Number_from_Action(testing_types, testing_number, action);

            this->Results.testing_types[t].Add(double(testing_types));
            this->Results.testing_number[t].Add(double(testing_number));

            this->Results.all_testing_types.push_back(testing_types);
            this->Results.all_testing_numbers.push_back(testing_number);

            total_testing_number += testing_number;

            if (this->ExpParams.debug_mode) {
                double numerator = 0.0, denominator = 0.0;
                Real.Get_Rewards_Numerator_Denominator(numerator, denominator, *state, action);
                if (numerator != 0 || denominator != 0) {
                    Results.Reward_Numerator.Add(numerator);
                    Results.Reward_Denomiantor.Add(denominator);
                }

            }
        }


//        // Debug code:
        if (isnan(reward)) abort();

        if (SearchParams.Verbose >= 1)
        {
            Real.DisplayAction(action, cout);
            Real.DisplayState(*state, cout);
            Real.DisplayObservation(*state, observation, cout);
            Real.DisplayReward(reward, cout);
        }

        if (termination)
        {
            cout << "Terminated" << endl;
            break;
        }
        outOfParticles = !mcts.Update(action, observation, reward);
        if (outOfParticles)
            break;

        if (timer.elapsed() > ExpParams.TimeOut)
        {
            cout << "Timed out after " << t << " steps in "
                 << Results.Time.GetTotal() << "seconds" << endl;
            break;
        }
    }

    if (outOfParticles)
    {
        cout << "Out of particles, finishing episode with SelectRandom" << endl;
        HISTORY history = mcts.GetHistory();
        while (++t < ExpParams.NumSteps)
        {
            int observation;
            double reward;

            // This passes real state into simulator!
            // SelectRandom must only use fully observable state
            // to avoid "cheating"
            int action = Simulator.SelectRandom(*state, history, mcts.GetStatus());
            termination = Real.Step(*state, action, observation, reward, cout);

            Results.Reward.Add(reward);
            undiscountedReturn += reward;
            discountedReturn += reward * discount;
            discount *= Real.GetDiscount();

            // Self implemented code:
            if (this->ExpParams.problems_name == "covid") {
                POMCP_STATE& pomcp_state = safe_cast<POMCP_STATE&>(*state);
                this->Results.S[t].Add(pomcp_state.S);
                this->Results.E[t].Add(pomcp_state.E);
                this->Results.I1[t].Add(pomcp_state.I1);
                this->Results.I2[t].Add(pomcp_state.I2);
                this->Results.H[t].Add(pomcp_state.H);
                this->Results.R[t].Add(pomcp_state.R);

                this->Results.all_S.push_back(pomcp_state.S);
                this->Results.all_E.push_back(pomcp_state.E);
                this->Results.all_I1.push_back(pomcp_state.I1);
                this->Results.all_I2.push_back(pomcp_state.I2);
                this->Results.all_H.push_back(pomcp_state.H);
                this->Results.all_R.push_back(pomcp_state.R);

                double testing_types = 0;
                int testing_number = 0;

                Real.Get_testing_Number_from_Action(testing_types, testing_number, action);
                this->Results.testing_types[t].Add(double(testing_types));
                this->Results.testing_number[t].Add(double(testing_number));

                this->Results.all_testing_types.push_back(testing_types);
                this->Results.all_testing_numbers.push_back(testing_number);

                total_testing_number += testing_number;

                if (this->ExpParams.debug_mode) {
                    double numerator = 0.0, denominator = 0.0;
                    Real.Get_Rewards_Numerator_Denominator(numerator, denominator, *state, action);
                    if (numerator != 0 || denominator != 0) {
                        Results.Reward_Numerator.Add(numerator);
                        Results.Reward_Denomiantor.Add(denominator);
                    }
                }
            }

            if (SearchParams.Verbose >= 1)
            {
                Real.DisplayAction(action, cout);
                Real.DisplayState(*state, cout);
                Real.DisplayObservation(*state, observation, cout);
                Real.DisplayReward(reward, cout);
            }

            if (termination)
            {
                cout << "Terminated" << endl;
                break;
            }

            history.Add(action, observation);
        }
    }


    Results.Time.Add(timer.elapsed());
    Results.UndiscountedReturn.Add(undiscountedReturn);
    Results.DiscountedReturn.Add(discountedReturn);


    double cum_I;
    if (this->ExpParams.problems_name == "covid") {
        const POMCP_STATE &pomcp_state = safe_cast<const POMCP_STATE &>(*state);
        cum_I = pomcp_state.I1 + pomcp_state.I2 + pomcp_state.R + pomcp_state.H;
        Results.Cumulative_I.Add(cum_I);
        Results.Total_Testing_Number.Add(total_testing_number);
    }

    cout << "Discounted return = " << discountedReturn
         << ", average = " << Results.DiscountedReturn.GetMean() << endl;
    cout << "Undiscounted return = " << undiscountedReturn
         << ", average = " << Results.UndiscountedReturn.GetMean() << endl;

    if (this->ExpParams.problems_name == "covid"){
        cout << "Cumulative_I = " << cum_I
             << ", average = " << Results.Cumulative_I.GetMean() << endl;
    }


}

void EXPERIMENT::MultiRun()
{
//#pragma omp parallel for
    for (int n = 0; n < ExpParams.NumRuns; n++)
    {
        cout << "Starting run " << n + 1 << " with "
             << SearchParams.NumSimulations << " simulations... " << endl;

        Run();

        if (Results.Time.GetTotal() > ExpParams.TimeOut)
        {
            cout << "Timed out after " << n << " runs in "
                 << Results.Time.GetTotal() << "seconds" << endl;
            break;
        }
    }
}

void EXPERIMENT::DiscountedReturn()
{
    cout << "Main runs" << endl;
    if (this->ExpParams.problems_name == "covid"){
        OutputFile << "Simulations,Runs,Undiscounted_return,Undiscounted_error,Discounted_return,Discounted_error,Time,Cumulative_I,Testing_Number_Total,"
                      "Testing_Number_Min,Testing_Number_Max,Testing_Number_Variance";

        processOutputFile << "Step,S,E,I1,I2,H,R,Testing_Types,Testing_Number" << endl;

        allOutputFile << "Run,Step,S,E,I1,I2,H,R,Testing_Types,Testing_Number" << endl;

        if (this->ExpParams.debug_mode){
            OutputFile << ",Reward_Numerator,Reward_Denominator"<< endl;
        } else {
            OutputFile << endl;
        }
    }else{
        OutputFile << "Simulations,Runs,Undiscounted_return,Undiscounted_error,Discounted_return,Discounted_error,Time\n";
    }

    SearchParams.MaxDepth = Simulator.GetHorizon(ExpParams.Accuracy, ExpParams.UndiscountedHorizon);
    ExpParams.SimSteps = Simulator.GetHorizon(ExpParams.Accuracy, ExpParams.UndiscountedHorizon);
    ExpParams.NumSteps = Real.GetHorizon(ExpParams.Accuracy, ExpParams.UndiscountedHorizon);

    for (int i = ExpParams.MinDoubles; i <= ExpParams.MaxDoubles; i++)
    {
        SearchParams.NumSimulations = 1 << i;
        SearchParams.NumStartStates = 1 << i;
        if (i + ExpParams.TransformDoubles >= 0)
            SearchParams.NumTransforms = 1 << (i + ExpParams.TransformDoubles);
        else
            SearchParams.NumTransforms = 1;
        SearchParams.MaxAttempts = SearchParams.NumTransforms * ExpParams.TransformAttempts;

        Results.Clear(); // Original Code. Clear the results before applying the new MultiRun function.
        //// User added code. For covid algorithm.
        if (this->ExpParams.problems_name=="covid") {
            Results.S = new STATISTIC[ExpParams.NumSteps];
            Results.E = new STATISTIC[ExpParams.NumSteps];
            Results.I1 = new STATISTIC[ExpParams.NumSteps];
            Results.I2 = new STATISTIC[ExpParams.NumSteps];
            Results.H = new STATISTIC[ExpParams.NumSteps];
            Results.R = new STATISTIC[ExpParams.NumSteps];
            Results.testing_types = new STATISTIC[ExpParams.NumSteps];
            Results.testing_number = new STATISTIC[ExpParams.NumSteps];
        }
        MultiRun();

        if (this->ExpParams.problems_name=="covid"){
            cout << "Simulations = " << SearchParams.NumSimulations << endl
                 << "Runs = " << Results.Time.GetCount() << endl
                 << "Undiscounted return = " << Results.UndiscountedReturn.GetMean()
                 << " +- " << Results.UndiscountedReturn.GetStdErr() << endl
                 << "Discounted return = " << Results.DiscountedReturn.GetMean()
                 << " +- " << Results.DiscountedReturn.GetStdErr() << endl
                 << "Cumulative I = " << Results.Cumulative_I.GetMean() << endl
                 << "Total_Testing_Number = " << Results.Total_Testing_Number.GetMean() << endl
                 << "Time = " << Results.Time.GetMean() << endl;
            OutputFile << SearchParams.NumSimulations << ", "
                       << Results.Time.GetCount() << ", "
                       << Results.UndiscountedReturn.GetMean() << ", "
                       << Results.UndiscountedReturn.GetStdErr() << ", "
                       << Results.DiscountedReturn.GetMean() << ", "
                       << Results.DiscountedReturn.GetStdErr() << ", "
                       << Results.Time.GetMean() << ", "
                       << Results.Cumulative_I.GetMean() << ", "
                       << Results.Total_Testing_Number.GetMean() << ", "
                       << Results.Total_Testing_Number.GetMin() << ", "
                       << Results.Total_Testing_Number.GetMax() << ", "
                       << Results.Total_Testing_Number.GetVariance();

            for (int debug_output_iter = 0; debug_output_iter < ExpParams.NumSteps; debug_output_iter++){
                processOutputFile << debug_output_iter << ", " << Results.S[debug_output_iter].GetMean() << ", " << Results.E[debug_output_iter].GetMean() << ", "
                                  << Results.I1[debug_output_iter].GetMean() << ", " << Results.I2[debug_output_iter].GetMean() << ", "
                                  << Results.H[debug_output_iter].GetMean()  << ", " << Results.R[debug_output_iter].GetMean() << ", " << Results.testing_types[debug_output_iter].GetMean() << ", "
                                  << Results.testing_number[debug_output_iter].GetMean()
                                  << endl;
            }

            for (int run_index = 0; run_index < ExpParams.NumRuns; run_index++ ){

                allOutputFile << run_index << ", " << 0 << ", " << Results.initial_cum_S
                              << ", " << Results.initial_cum_E << ", " << Results.initial_cum_I1 << ", " << Results.initial_cum_I2
                              << ", " << Results.initial_cum_H << ", " << Results.initial_cum_R
                              << ", " << 1 << ", " << 0
                              << endl;

                for (int step_index = 0; step_index < ExpParams.NumSteps; step_index++){
                    int current_index = ExpParams.NumSteps * (run_index) + step_index;
                    allOutputFile << run_index << ", " << step_index+1 << ", " << Results.all_S[current_index]
                    << ", " << Results.all_E[current_index] << ", " << Results.all_I1[current_index] << ", " << Results.all_I2[current_index]
                    << ", " << Results.all_H[current_index] << ", " << Results.all_R[current_index]
                    << ", " << Results.all_testing_types[current_index] << ", " << Results.all_testing_numbers[current_index]
                    << endl;
                }
            }

            if (this->ExpParams.debug_mode){
                OutputFile << ", " << Results.Reward_Numerator.GetMean() << ", "
                           << Results.Reward_Denomiantor.GetMean() << endl;
            } else {
                OutputFile << endl;
            }
        } else {
            cout << "Simulations = " << SearchParams.NumSimulations << endl
                 << "Runs = " << Results.Time.GetCount() << endl
                 << "Undiscounted return = " << Results.UndiscountedReturn.GetMean()
                 << " +- " << Results.UndiscountedReturn.GetStdErr() << endl
                 << "Discounted return = " << Results.DiscountedReturn.GetMean()
                 << " +- " << Results.DiscountedReturn.GetStdErr() << endl
                 << "Time = " << Results.Time.GetMean() << endl;
            OutputFile << SearchParams.NumSimulations << ", "
                       << Results.Time.GetCount() << ", "
                       << Results.UndiscountedReturn.GetMean() << ", "
                       << Results.UndiscountedReturn.GetStdErr() << ", "
                       << Results.DiscountedReturn.GetMean() << ", "
                       << Results.DiscountedReturn.GetStdErr() << ", "
                       << Results.Time.GetMean() << endl;
        }
    }
}

void EXPERIMENT::AverageReward()
{
    cout << "Main runs" << endl;
    if (this->ExpParams.problems_name == "covid"){
        OutputFile << "Simulations,Steps,Average_reward,Average_time,Cumulative_I,Testing_Number_Total,"
                      "Testing_Number_Min,Testing_Number_Max,Testing_Number_Variance";
        if (this->ExpParams.debug_mode){
            OutputFile << ",Reward_Numerator,Reward_Denominator"<< endl;

            processOutputFile << "Step,S,E,I1,I2,H,R,Testing_Types,Testing_Number" << endl;

            allOutputFile << "Run,Step,S,E,I1,I2,H,R,Testing_Types,Testing_Number" << endl;
        } else {
            OutputFile << endl;
        }
    }else {
        OutputFile << "Simulations,Steps,Average_reward,Average_time"<<endl;
    }

    SearchParams.MaxDepth = Simulator.GetHorizon(ExpParams.Accuracy, ExpParams.UndiscountedHorizon);
    ExpParams.SimSteps = Simulator.GetHorizon(ExpParams.Accuracy, ExpParams.UndiscountedHorizon);

    for (int i = ExpParams.MinDoubles; i <= ExpParams.MaxDoubles; i++)
    {
        SearchParams.NumSimulations = 1 << i;
        SearchParams.NumStartStates = 1 << i;
        if (i + ExpParams.TransformDoubles >= 0)
            SearchParams.NumTransforms = 1 << (i + ExpParams.TransformDoubles);
        else
            SearchParams.NumTransforms = 1;
        SearchParams.MaxAttempts = SearchParams.NumTransforms * ExpParams.TransformAttempts;

        Results.Clear();
        if (this->ExpParams.problems_name=="covid") {
            Results.S = new STATISTIC[ExpParams.NumSteps];
            Results.E = new STATISTIC[ExpParams.NumSteps];
            Results.I1 = new STATISTIC[ExpParams.NumSteps];
            Results.I2 = new STATISTIC[ExpParams.NumSteps];
            Results.H = new STATISTIC[ExpParams.NumSteps];
            Results.R = new STATISTIC[ExpParams.NumSteps];
            Results.testing_number = new STATISTIC[ExpParams.NumSteps];
        }
        Run();

        if (this->ExpParams.problems_name == "covid"){
            cout << "Simulations = " << SearchParams.NumSimulations << endl
                 << "Steps = " << Results.Reward.GetCount() << endl
                 << "Average reward = " << Results.Reward.GetMean()
                 << " +- " << Results.Reward.GetStdErr() << endl
                 << "Average time = " << Results.Time.GetMean() / Results.Reward.GetCount() << endl
                 << "Avg_Testing_Number = " << Results.Total_Testing_Number.GetMean() << endl
                 << "Cummulative_I = " << Results.Cumulative_I.GetMean() << endl;
            OutputFile << SearchParams.NumSimulations << ", "
                       << Results.Reward.GetCount() << ", "
                       << Results.Reward.GetMean() << ", "
                       << Results.Reward.GetStdErr() << ", "
                       << Results.Time.GetMean() / Results.Reward.GetCount() << ","
                       << Results.Cumulative_I.GetMean() << ","
                       << Results.Total_Testing_Number.GetMean() << ", "
                       << Results.Total_Testing_Number.GetMin() << ", "
                       << Results.Total_Testing_Number.GetMax() << ", "
                       << Results.Total_Testing_Number.GetVariance();

            for (int debug_output_iter = 0; debug_output_iter < ExpParams.NumSteps; debug_output_iter++){
                processOutputFile << debug_output_iter << ", " << Results.S[debug_output_iter].GetMean() << ", " << Results.E[debug_output_iter].GetMean() << ", "
                                  << Results.I1[debug_output_iter].GetMean() << ", " << Results.I2[debug_output_iter].GetMean() << ", "
                                  << Results.H[debug_output_iter].GetMean()  << ", " << Results.R[debug_output_iter].GetMean() << ", " << Results.testing_types[debug_output_iter].GetMean() << ", "
                                  << Results.testing_number[debug_output_iter].GetMean()
                                  << endl;
            }

            for (int run_index = 0; run_index < ExpParams.NumRuns; run_index++ ){

                allOutputFile << run_index << ", " << 0 << ", " << Results.initial_cum_S
                              << ", " << Results.initial_cum_E << ", " << Results.initial_cum_I1 << ", " << Results.initial_cum_I2
                              << ", " << Results.initial_cum_H << ", " << Results.initial_cum_R
                              << ", " << 1 << ", " << 0
                              << endl;

                for (int step_index = 0; step_index < ExpParams.NumSteps; step_index++){
                    int current_index = ExpParams.NumSteps * (run_index) + step_index;
                    allOutputFile << run_index << ", " << step_index+1 << ", " << Results.all_S[current_index]
                                  << ", " << Results.all_E[current_index] << ", " << Results.all_I1[current_index] << ", " << Results.all_I2[current_index]
                                  << ", " << Results.all_H[current_index] << ", " << Results.all_R[current_index]
                                  << ", " << Results.all_testing_types[current_index] << ", " << Results.all_testing_numbers[current_index]
                                  << endl;
                }
            }

            if (this->ExpParams.debug_mode){
                OutputFile << ", " << Results.Reward_Numerator.GetMax() << ", "
                           << Results.Reward_Denomiantor.GetMax()<< endl;
            } else {
                OutputFile << endl;
            }
        } else {
            cout << "Simulations = " << SearchParams.NumSimulations << endl
                 << "Steps = " << Results.Reward.GetCount() << endl
                 << "Average reward = " << Results.Reward.GetMean()
                 << " +- " << Results.Reward.GetStdErr() << endl
                 << "Average time = " << Results.Time.GetMean() / Results.Reward.GetCount() << endl;
            OutputFile << SearchParams.NumSimulations << ", "
                       << Results.Reward.GetCount() << ", "
                       << Results.Reward.GetMean() << ", "
                       << Results.Reward.GetStdErr() << ", "
                       << Results.Time.GetMean() / Results.Reward.GetCount() << endl;
        }
    }
}

//----------------------------------------------------------------------------
