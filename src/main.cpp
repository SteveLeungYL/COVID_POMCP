//#include "game.h"
#include "pomcp.h"
#include "battleship.h"
#include "mcts.h"
#include "network.h"
#include "pocman.h"
#include "rocksample.h"
#include "tag.h"
#include "experiment.h"
#include <boost/program_options.hpp>
#include<fstream>
#include<string>

using namespace std;
using namespace boost::program_options;

void UnitTests()
{
    cout << "Testing UTILS" << endl;
    UTILS::UnitTest();
    cout << "Testing COORD" << endl;
    COORD::UnitTest();
    cout << "Testing MCTS" << endl;
    MCTS::UnitTest();
}

void disableBufferedIO(void)
{
    setbuf(stdout, NULL);
    setbuf(stdin, NULL);
    setbuf(stderr, NULL);
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stdin, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
}

vector<vector<int>*>* enumeratePuredStrategies(int maximum_testing_number, int testing_group_number, int testing_action_step)
{
    vector<vector<int>*> *testing_strategies = new vector<vector<int>*>();

    for (int asymptomatic_testing_number = 0; asymptomatic_testing_number <= maximum_testing_number; asymptomatic_testing_number+=testing_action_step){
        vector<int>* testing_strategy = new vector<int>();

        testing_strategy->push_back(asymptomatic_testing_number);
        int symptomatic_testing_number = maximum_testing_number - asymptomatic_testing_number;
        testing_strategy->push_back(symptomatic_testing_number);
        testing_strategy->push_back(testing_group_number);

        testing_strategies->push_back(testing_strategy);
    }

    return testing_strategies;
}


int main(int argc, char* argv[])
{
    MCTS::PARAMS searchParams;
    EXPERIMENT::PARAMS expParams;
    SIMULATOR::KNOWLEDGE knowledge;
    string problem, outputfile, policy;
    int size, number, treeknowledge = 1, rolloutknowledge = 1, smarttreecount = 10;
    double smarttreevalue = 1.0;

    //user defind vars. Used for COVID-19 testing policy generator.
    int maximum_testing_number = 200;
    int testing_group_number = 32;
    int testing_action_step = 10;

    map<int, vector<int>*> action_ID_to_action_map;
    map<vector<int>*, int> action_to_action_ID_map;

    int total_population = 50000;
    double R0 = 2.0;
    int time_step = 1;
//    int total_time_step = 100;
    double testing_sensitivity = 0.90, testing_specitivity = 0.998, proportion_of_I2_get_symptomatic_testing = 0.3, asymptomatic_rate = 0.7;
    double I_class_ratio = 0.005;
    double R_class_ratio = 0.005;
    double symptomatic_testing_COVID_rate = 0.76;
    int fixed_action_ID = 0;
    int is_random_action = 0;
    int debug_mode = 0;

    options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("test", "run unit tests")

            //////// User defined inputs.
            ("total_population", value<int>(&total_population), "determine the total population of the target testing region.")
            ("maximum_testing_number", value<int>(&maximum_testing_number), "determine the maximum testing methods, include both Asymptomatic testing and Symptomatic testing.")
//            ("maximum_testing_group_number", value<int>(&maximum_testing_group_number), "determine the maximum testing group number, include only for Asymptomatic testing.")
            ("testing_group_number", value<int>(&testing_group_number), "determine the current testing group number. works only for Asymptomatic testing.")
            ("testing_action_step", value<int>(&testing_action_step), "determine the step of the action range. Action = {0, maximum_testing_number, action_step}, will impact the action number. Default 1.")
            ("R0_value", value<double>(&R0),"determine the R0 value in the target testing region (Default 2.0)")
            ("time_step", value<int>(&time_step), "determine how many time steps between each action changes (Default 1)")
//            ("total_time_step", value<int>(&total_time_step), "determine the total time steps that the simulation SEIR model will run (Default 100)")
            ("testing_sensitivity", value<double>(&testing_sensitivity), "determine the testing sensitivity, include asymptomatic testing and symptomatic"
                                                                         "testing (Default 0.90)")
            ("testing_specitivity", value<double>(&testing_specitivity), "determine the testing specitivity, include asymptomatic testing and symptomatic"
                                                                         "testing (Default 0.998)")
            ("proportion_of_I2_get_symptomatic_testing", value<double>(&proportion_of_I2_get_symptomatic_testing), "determine the proportion of patients that are "
                                                                                                                   "symptomatic get symptomatic_testing (Default 0.5)")
            ("fix_action_ID", value<int>(&fixed_action_ID), "fix the action that the POMCP model would take. (For Debug Purpose)")
            ("is_random_action", value<int>(&is_random_action), "determine whether we use random action.")
            ("debug_mode", value<int>(&debug_mode), "determine whether we enter debug mode. (For POMCP problem only) (Output Reward_Numerator and Reward_Denominator)")
            ("I_class_ratio", value<double>(&I_class_ratio), "determine the ratio of the I class over the total population in the start state.")
            ("R_class_ratio", value<double>(&R_class_ratio), "determine the ratio of the R class over the total population in the start state.")
            ("asymptomatic_rate", value<double>(&asymptomatic_rate), "determine the ratio of the asymptomatic cases in the total infectious class.")
            ("symptomatic_testing_COVID_rate", value<double>(&symptomatic_testing_COVID_rate), "determine the positive rate of symptomatic testing.")
            /////// Dr. Yadav's inputs.
//    ("gamefile", value<string>(&gamefile), "file containing game payoffs")
//    ("MSD", value<double>(&K), "discretization of the finite mixed strategy space")
//    ("LambdaDis", value<double>(&LambdaDis), "discretization of the lambda space")
//    ("maxLambda", value<double>(&maxLambda), "max lambda value")
            ("problem", value<string>(&problem), "problem to run")
            ("outputfile", value<string>(&outputfile)->default_value("output.csv"), "summary output file")
            ("policy", value<string>(&policy), "policy file (explicit POMDPs only)")
            ("size", value<int>(&size), "size of problem (problem specific)")
            ("number", value<int>(&number), "number of elements in problem (problem specific)")
            ("timeout", value<double>(&expParams.TimeOut), "timeout (seconds)")
            ("mindoubles", value<int>(&expParams.MinDoubles), "minimum power of two simulations")
            ("maxdoubles", value<int>(&expParams.MaxDoubles), "maximum power of two simulations")
            ("runs", value<int>(&expParams.NumRuns), "number of runs")
            ("accuracy", value<double>(&expParams.Accuracy), "accuracy level used to determine horizon")
            ("horizon", value<int>(&expParams.UndiscountedHorizon), "horizon to use when not discounting")
            ("num_steps", value<int>(&expParams.NumSteps), "number of steps to run when using average reward")
            ("verbose", value<int>(&searchParams.Verbose), "verbosity level")
            ("autoexploration", value<bool>(&expParams.AutoExploration), "Automatically assign UCB exploration constant")
            ("exploration", value<double>(&searchParams.ExplorationConstant), "Manual value for UCB exploration constant")
            ("usetransforms", value<bool>(&searchParams.UseTransforms), "Use transforms")
            ("transformdoubles", value<int>(&expParams.TransformDoubles), "Relative power of two for transforms compared to simulations")
            ("transformattempts", value<int>(&expParams.TransformAttempts), "Number of attempts for each transform")
            ("userave", value<bool>(&searchParams.UseRave), "RAVE")
            ("ravediscount", value<double>(&searchParams.RaveDiscount), "RAVE discount factor")
            ("raveconstant", value<double>(&searchParams.RaveConstant), "RAVE bias constant")
            ("treeknowledge", value<int>(&knowledge.TreeLevel), "Knowledge level in tree (0=Pure, 1=Legal, 2=Smart)")
            ("rolloutknowledge", value<int>(&knowledge.RolloutLevel), "Knowledge level in rollouts (0=Pure, 1=Legal, 2=Smart)")
            ("smarttreecount", value<int>(&knowledge.SmartTreeCount), "Prior count for preferred actions during smart tree search")
            ("smarttreevalue", value<double>(&knowledge.SmartTreeValue), "Prior value for preferred actions during smart tree search")
            ("disabletree", value<bool>(&searchParams.DisableTree), "Use 1-ply rollout action selection")
            ;

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    expParams.problems_name = problem;
    expParams.fixed_action_ID = fixed_action_ID;
    expParams.is_random_action = is_random_action;
    expParams.debug_mode = debug_mode;

    if (vm.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("problem") == 0)
    {
        cout << "No problem specified" << endl;
        return 1;
    }

    if (vm.count("test"))
    {
        cout << "Running unit tests" << endl;
        UnitTests();
        return 0;
    }

    SIMULATOR* real = 0;
    SIMULATOR* simulator = 0;

    if (problem == "battleship")
    {
        real = new BATTLESHIP(size, size, number);
        simulator = new BATTLESHIP(size, size, number);
    }
    else if (problem == "pocman")
    {
        real = new FULL_POCMAN;
        simulator = new FULL_POCMAN;
    }
    else if (problem == "network")
    {
        real = new NETWORK(size, number);
        simulator = new NETWORK(size, number);
    }
    else if (problem == "rocksample")
    {
        real = new ROCKSAMPLE(size, number);
        simulator = new ROCKSAMPLE(size, number);
    }
    else if (problem == "tag")
    {
        real = new TAG(number);
        simulator = new TAG(number);
    }

    else if (problem=="covid"){
        assert(total_population);

        vector<vector<int>*>* all_pured_strategies = enumeratePuredStrategies(maximum_testing_number, testing_group_number, testing_action_step);
        for (int i = 0; i < all_pured_strategies->size(); i++){
            action_ID_to_action_map[i] = all_pured_strategies->at(i);
            action_to_action_ID_map[all_pured_strategies->at(i)] = i;
        }
        // Both maps created.
        real = new POMCP(&action_ID_to_action_map, &action_to_action_ID_map, total_population,
                         R0, time_step, testing_sensitivity, testing_specitivity,
                         proportion_of_I2_get_symptomatic_testing, I_class_ratio, R_class_ratio, asymptomatic_rate, symptomatic_testing_COVID_rate);
        simulator = new POMCP(&action_ID_to_action_map, &action_to_action_ID_map, total_population,
                              R0, time_step, testing_sensitivity, testing_specitivity,
                              proportion_of_I2_get_symptomatic_testing, I_class_ratio, R_class_ratio, asymptomatic_rate, symptomatic_testing_COVID_rate);
    }
    else
    {
        cout << "Unknown problem" << endl;
        exit(1);
    }


    simulator->SetKnowledge(knowledge);
    EXPERIMENT experiment(*real, *simulator, outputfile, expParams, searchParams);
    experiment.DiscountedReturn();

    delete real;
    delete simulator;
    return 0;
}
