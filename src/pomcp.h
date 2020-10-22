//
//  POMCP.hpp
//  POMCP
//
//  Created by Yu Liang on 6/25/20.
//  Copyright © 2020 Yu Liang. All rights reserved.
//

#ifndef POMCP_hpp
#define POMCP_hpp

#include <stdio.h>

#include "simulator.h"
#include<vector>
#include<map>


using namespace std;

class POMCP_STATE : public STATE
{
public:
    double S, E, I1, I2, H, R;
//    int T;
//    int total_testing_number;
};

class POMCP : public SIMULATOR
{
public:
    ///constructor
    POMCP(map<int, vector<int> *> *action_ID_to_action_map, map<vector<int> *, int> *action_to_action_ID_map,
          int total_population, double R0, int time_step, double testing_sensitivity, double testing_specificity,
          double proportion_of_I2_get_symptomatic_testing, double I_class_ratio, double R_class_ratio,
          double asymptomatic_rate, double symptomatic_testing_COVID_rate);

    virtual STATE* CreateStartState() const;
    
    // Free memory for state
    virtual void FreeState(STATE* state) const;
    
    // Update state according to action, and get observation and reward.
    // Return value of true indicates termination of episode (if episodic)
    virtual bool Step(STATE& state, int action,
                      int& observation, double& reward, std::ostream& ostr) const;
    
    // Create new state and copy argument (must be same type)
    virtual STATE* Copy(const STATE& state) const;
    
    // Sanity check
    virtual void Validate(const STATE& state) const;
    
    // Modify state stochastically to some related state
//    virtual bool LocalMove(STATE& state, const HISTORY& history,
//                           int stepObs, const STATUS& status) const;
    
    
    //    int SelectRandom(const STATE& state, const HISTORY& history,
    //        const STATUS& status) const;
    
    // Generate set of legal actions
    virtual void GenerateLegal(const STATE& state, const HISTORY& history,
                               std::vector<int>& actions, const STATUS& status) const;
    
    virtual void DisplayBeliefs(const BELIEF_STATE& beliefState,
                                std::ostream& ostr) const;
    virtual void DisplayState(const STATE& state, std::ostream& ostr) const;
    virtual void DisplayAction(int action, std::ostream& ostr) const;
    virtual void DisplayObservation(const STATE& state, int observation, std::ostream& ostr) const;

    //// Self added virtual function.
    virtual void Get_testing_Number_from_Action(double& testing_type, int& testing_number, int action_ID) const;
    virtual int Get_Rewards_Numerator_Denominator(double &numerator, double &denominator, STATE &state, int action) const;

    void update_SEIR(double &S, double &E, double &I1, double &I2, double &R, double &H) const;

    void
    asymptomatic_testing(int testing_number, int testing_group_number, double &S, double &E, double &I1, double &I2,
                         double &R, double &H, int &observation) const;
    void symptomatic_testing(int testing_number, double &I2, double &H, int &observation) const;

private:
    mutable MEMORY_POOL<POMCP_STATE> MemoryPool;
    ////user defined variables
    
    map<int, vector<int>*>* action_ID_to_action_map;
    map<vector<int>*, int>* action_to_action_ID_map;
    int total_population = 0;
//    int maximum_testing_number = 0;
//    int maximum_testing_group_number = 0;
    int time_step = 1;
//    int total_time_step = 0;
    double R0 = 1.0;

    // SEIR updates related:
    double beta;
    double theta;
    double gamma_I1_R; // 14 days infectious period.
    double gamma_I2_R; // 14 days infectious period.

    // Testing related
    double testing_sensitivity = 0.90;
    double testing_specificity = 0.90;
    double proportion_of_I2_get_symptomatic_testing = 0.2;
    double asymptomatic_rate = 0.7;
    double hospital_recovery_rate = 1.0 / 14.0;
    double symptomatic_testing_COVID_rate = 0.76;

    double I_class_ratio = 0.0, R_class_ratio = 0.0;
};
#endif
