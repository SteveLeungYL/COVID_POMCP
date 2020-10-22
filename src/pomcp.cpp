//
//  POMCP.cpp
//  POMCP
//
//  Created by Yu Liang on 6/25/20.
//  Copyright © 2020 Yu Liang. All rights reserved.
//

#include "pomcp.h"

#include<time.h>
//#include<math.h>
#include<string>
#include<iostream>
#include<fstream>
#include <cmath>
#include <random>

using namespace std;

POMCP::POMCP(map<int, vector<int> *> *action_ID_to_action_map, map<vector<int> *, int> *action_to_action_ID_map,
             int total_population, double R0, int time_step, double testing_sensitivity, double testing_specificity,
             double proportion_of_I2_get_symptomatic_testing, double I_class_ratio, double R_class_ratio,
             double asymptomatic_rate, double symptomatic_testing_COVID_rate)
{

    this->action_ID_to_action_map = action_ID_to_action_map;
    this->action_to_action_ID_map = action_to_action_ID_map;
    this->total_population = total_population;
//    this->maximum_testing_number = maximum_testing_number;
//    this->maximum_testing_group_number = maximum_testing_group_number;
    this->R0 = R0;
    this->time_step = time_step;
//    this->total_time_step = total_time_step;
    this->testing_sensitivity = testing_sensitivity;
    this->testing_specificity = testing_specificity;
    this->proportion_of_I2_get_symptomatic_testing = proportion_of_I2_get_symptomatic_testing;
    this->I_class_ratio = I_class_ratio;
    this->R_class_ratio = R_class_ratio;
    this->asymptomatic_rate = asymptomatic_rate;
    this->symptomatic_testing_COVID_rate = symptomatic_testing_COVID_rate;

    theta = 1.0/5.0;
    gamma_I1_R = 1.0/14.0; // 14 days infectious period.
    gamma_I2_R = 1.0/14.0; // 14 days infectious period.

    beta = R0 * (gamma_I1_R + gamma_I2_R) / 2.0;


    /////variables from superclass
    NumActions = action_ID_to_action_map->size();
    NumObservations = total_population;
    Discount = 1.0;
    RewardRange = 0.1;////redundant variable as long as autoexploration is false

}

STATE* POMCP::CreateStartState() const
{
    /* E = 6/14 * I. About 0.77% to 1.59% of population.
     * I = 1.8% ~ 3.7% of population.
     * R = 6.25% ~ 12.5% of population.
     */

//    assert(this->total_population);
//    assert(this->total_time_step);

    POMCP_STATE* state = MemoryPool.Allocate();

//    random_device rd;
//    mt19937_64 gen(rd());
//    uniform_real_distribution<float> rand_dis(0.0, 1.0);
//    double temp_random = rand_dis(gen);

    double I = (double)total_population * I_class_ratio;
    double I1 = I * this->asymptomatic_rate;
    double I2 = I - I1;
    double E = 6.0/14.0 * I;

//    temp_random = rand_dis(gen);
    double R = (double)total_population * R_class_ratio;

    double S = (double)total_population - E - I - R;

    state->S = S;
    state->E = E;
    state->I1 = I1;
    state->I2 = I2;

    state->H = 0;

    state->R = R;
//    state->T = this->total_time_step;
//    state->total_testing_number = 0;


//    cout<<"SEIRT: " << S << " " << E << " " << I1 << " " << I2 << " " << R << endl;

    return state;
}

void POMCP::FreeState(STATE *state) const
{
    POMCP_STATE* bsstate = safe_cast<POMCP_STATE*>(state);
    MemoryPool.Free(bsstate);
}

STATE* POMCP::Copy(const STATE& state) const
{
    assert(state.IsAllocated());
    const POMCP_STATE& oldstate = safe_cast<const POMCP_STATE&>(state);
    POMCP_STATE* newstate = MemoryPool.Allocate();
    *newstate = oldstate;
    return newstate;
}

void POMCP::Validate(const STATE& state) const  // how to varified the isValid status?
{
    bool isValid = true;
    const POMCP_STATE& pomcp_state = safe_cast<const POMCP_STATE&>(state);

    double S, E, I1, I2, H, R;
//    int T, current_total_testing_number;

    S = pomcp_state.S;
    E = pomcp_state.E;
    I1 = pomcp_state.I1;
    I2 = pomcp_state.I2;
    H = pomcp_state.H;
    R = pomcp_state.R;
//    T = pomcp_state.T;
//    current_total_testing_number = pomcp_state.total_testing_number;

    if (S > total_population || E > total_population || I1 > total_population || I2 > total_population || R > total_population || H > total_population) isValid = false;
    else isValid = true;

    assert(isValid);
}

void POMCP::DisplayBeliefs(const BELIEF_STATE& beliefState,
                           ostream& ostr) const
{
    const POMCP_STATE *pomcp_belief_state;

    ostr<<"Belief: \n";
    for (int i=0;i<beliefState.GetNumSamples();i++)
    {
        pomcp_belief_state = safe_cast<const POMCP_STATE*>(beliefState.GetSample(i));
        ostr << "Belief State: S:" << pomcp_belief_state->S << "; E: " << pomcp_belief_state->E << "; I1: "
             << pomcp_belief_state->I1 << "; I2: " << pomcp_belief_state->I2<< "; R: " << pomcp_belief_state->R
             << "; H: " << pomcp_belief_state->H << ". \n";
    }
}

void POMCP::DisplayState(const STATE& state, ostream& ostr) const
{
    const POMCP_STATE& pomcp_state = safe_cast<const POMCP_STATE&>(state);
    ostr<<"State: S:"<<pomcp_state.S<<"; E: "<< pomcp_state.E <<"; I1: " <<pomcp_state.I1 <<"; I2: "<<pomcp_state.I2
        << "; R: "<<pomcp_state.R << "; H: "<< pomcp_state.H << ". " << endl;
}

void POMCP::DisplayObservation(const STATE& state, int observation, ostream& ostr) const
{
    ostr<<"Observation: Current Testing Positive number: "<<observation<<"\n";
}

void POMCP::DisplayAction(int action, ostream& ostr) const
{
    ostr<<"Action: \n";
    vector<int>* action_details = action_ID_to_action_map->find(action)->second;

    if (action_details->at(0) == 0 && action_details->size() == 3 && action_details->at(2) != 0){
        ostr<<"Action: Asymptomatic Testing. \n";
        ostr<<"Testing number per time step: "<<action_details->at(1) <<". Testing number per testing group: " << action_details->at(2) <<". \n";
    }
    else if(action_details->at(0) == 1 && action_details->size() == 2){
        ostr<<"Action: Symptomatic Testing. \n";
        ostr<<"Testing number per time step: "<<action_details->at(1)<<". \n";
    } else {
        cerr << "Invalid action produced. Current Action ID: " << action << ". Exit!";
        assert(false);
    }
    return;
}

//bool POMCP::LocalMove(STATE& state, const HISTORY& history, int stepObs, const STATUS& status) const
//{
//    ///we do nothing assuming that a state is stoschatically related to another state if the 2nd state can be reached from 1st by transition
//    //so we do not change state
//    ///however this may lead to particle deprivation
////    cout<<"Readding same particle\n";
////    return true;
//}

void POMCP::GenerateLegal(const STATE& state, const HISTORY& history, vector<int>& actions, const STATUS& status) const
{
    actions.clear();
    for (map<int, vector<int>*>::iterator it = this->action_ID_to_action_map->begin(); it != this->action_ID_to_action_map->end(); it++)
    {
        actions.push_back(it->first);
    }
}

void POMCP::update_SEIR(double &S, double &E, double &I1, double &I2, double &R, double &H) const {

    double ori_S, ori_E, ori_I1, ori_I2, ori_R, ori_H;
    ori_S = S;
    ori_E = E;
    ori_I1 = I1;
    ori_I2 = I2;
    ori_R = R;
    ori_H = H;

    double alpha_1 = this->asymptomatic_rate;
    double alpha_2 = 1.0 - alpha_1;
    double omega = hospital_recovery_rate;

    S = ori_S - (beta * (ori_I1 + ori_I2) * ori_S) / total_population;
    E = ori_E - theta * ori_E +  (beta * (ori_I1 + ori_I2) * ori_S) / total_population;
    I1 = ori_I1 - ori_I1 * gamma_I1_R + alpha_1 * theta * ori_E;
    I2 = ori_I2 - ori_I2 * gamma_I2_R + alpha_2 * theta * ori_E;

    H = ori_H - ori_H * omega;
    R = R + ori_I1 * gamma_I1_R + ori_I2 * gamma_I2_R + ori_H * omega;

}

void
POMCP::asymptomatic_testing(int testing_number, int testing_group_number, double &S, double &E, double &I1, double &I2,
                            double &R, double &H, int &observation) const{

//    double proportion_of_R = R / (S + E + I1 + R);

//    double ori_I1 = I1; // To calculate rewards.
//    double ori_I2 = I2; // To calculate rewards.

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> rand_dis(0.0, 1.0);

    for (int i = 0; i < testing_number; i++){


        double current_tested_I1 = 0;
        double current_tested_I2 = 0;

        double current_tested_I1_pos = 0; // True_positive.
        double current_tested_I2_pos = 0; // True_positive.

        if (S < 0.0) S = 0.0;
        if (E < 0.0) S = 0.0;
        if (I1 < 0.0) S = 0.0;
        if (I2 < 0.0) S = 0.0;
        if (R < 0.0) S = 0.0;
        if (H < 0.0) H = 0.0;

        double proportion_of_I1 = I1 / (S + E + I1 + I2 + R);
        double proportion_of_I2 = I2 / (S + E + I1 + I2 + R);

        for (int j = 0; j < testing_group_number; j++){

            double temp_random = (double)rand_dis(gen);

            if (temp_random < proportion_of_I1){
                double temp_random_2 = (double)rand_dis(gen);
                if (temp_random_2 <= testing_sensitivity) {
                    current_tested_I1 += 1.0;
                    current_tested_I1_pos += 1.0;
                }
                else current_tested_I1 += 1.0;
            }

            else if (temp_random < proportion_of_I1 + proportion_of_I2){
                double temp_random_2 = (double)rand_dis(gen);
                if (temp_random_2 <= testing_sensitivity) {
                    current_tested_I2+=1.0;
                    current_tested_I2_pos += 1.0;
                }
                else current_tested_I2 += 1.0;
            }
        }

        if(current_tested_I1_pos || current_tested_I2_pos){

            I1 -= current_tested_I1_pos;
            if (I1 <0) I1 = 0.0;

            I2 -= current_tested_I2_pos;
            if (I2 < 0) I2 = 0.0;

            H += current_tested_I1_pos + current_tested_I2_pos;
            observation += (int) (current_tested_I1_pos + current_tested_I2_pos);
        }
        testing_number -= (current_tested_I1_pos + current_tested_I2_pos) * log(testing_group_number) / log(2.0);
    }
}

void POMCP::symptomatic_testing(int testing_number, double &I2, double &H, int &observation) const{

    double ori_I2 = I2;
    int current_tested_I2_pos = 0;

    testing_number = int(double(testing_number) * this->symptomatic_testing_COVID_rate);

    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> rand_dis(0.0, 1.0);

    for (int i = 0; i < (int)ori_I2; i++){
        double temp_random = (double)rand_dis(gen);
        if ( temp_random <= (this->proportion_of_I2_get_symptomatic_testing * this->testing_sensitivity)){
            if (I2 > 0) {
                I2 -= 1;
                H += 1;
                current_tested_I2_pos += 1;
            } else{ // If I2 <= 0;
                I2 = 0;
                break;
            }
        }
        --testing_number;
        if (testing_number <= 0) break;
    }
    observation += current_tested_I2_pos;
    return;
}


bool POMCP::Step(STATE& state, int action, int& observation, double& reward, std::ostream& ostr) const
{
    POMCP_STATE& pomcp_state = safe_cast<POMCP_STATE&>(state);

    double S, E, I1, I2, R, H;
//    int T, current_total_testing_number;
    S = pomcp_state.S;
    E = pomcp_state.E;
    I1 = pomcp_state.I1;
    I2 = pomcp_state.I2;
    R = pomcp_state.R;
    H = pomcp_state.H;
//    current_total_testing_number = pomcp_state.total_testing_number;

//    cout<<"SEIRT: " << S << " " << E << " " << I1 << " " << I2 << " " << R << " " << T << " " << current_total_testing_number << endl;

    double ori_I = I1 + I2;
//    double I2_pos = 0.0;

    observation = 0;
    reward = 0.0;


    vector<int>* action_details = this->action_ID_to_action_map->at(action);
    int asymptomatic_testing_number = action_details->at(0);
    int symptomatic_testing_number = action_details->at(1);
    int testing_group_number = action_details->at(2);

    for (int iter = 0; iter < time_step; iter++){
        // Begin single step's SEIR updates.
        update_SEIR(S, E, I1, I2, R, H);

        asymptomatic_testing(asymptomatic_testing_number, testing_group_number, S, E, I1, I2, R, H,
                             observation);

        symptomatic_testing(symptomatic_testing_number, I2, H, observation);

    }



// Reward Function.
    reward = -(I1 + I2);
//    if (testing_types == 0){ // Asymptomatic testing.
//        reward /= max(double(testing_number), 0.99);
//    }
//    else{ //Symptomatic testing.
//        reward /= max(double(testing_number), 0.99);
//    }


    if (reward != reward) assert(false);
//    if (reward < 0.0) reward = 0.0;

    pomcp_state.S = S;
    pomcp_state.E = E;
    pomcp_state.I1 = I1;
    pomcp_state.I2 = I2;
    pomcp_state.R = R;
    pomcp_state.H = H;

    return false;  // Return false means algorithms not terminated!!!
}

void POMCP::Get_testing_Number_from_Action(double& testing_type, int& testing_number, int action_ID) const {

    vector<int>* action_details = this->action_ID_to_action_map->at(action_ID);
    testing_number = action_details->at(0);
    testing_number += action_details->at(1);

    testing_type = action_details->at(0) / (double)testing_number;
}

int POMCP::Get_Rewards_Numerator_Denominator(double &numerator, double &denominator, STATE &state, int action) const{

    // Setup required variables.
    POMCP_STATE& pomcp_state = safe_cast<POMCP_STATE&>(state);
    double S, E, I1, I2, R, H;
//    int T, current_total_testing_number;
    S = pomcp_state.S;
    E = pomcp_state.E;
    I1 = pomcp_state.I1;
    I2 = pomcp_state.I2;
    R = pomcp_state.R;
    H = pomcp_state.H;


//    vector<int>* action_details = this->action_ID_to_action_map->at(action);
//    int testing_types = action_details->at(0);
//    int testing_number = action_details->at(1);
//    int testing_group_number = 0;
//    if (testing_types == 0) {
//        testing_group_number = action_details->at(2);;
//    }


    numerator = S;
//    denominator = max(double(testing_number), 0.99);
    denominator = 0;

    return 0;
}

