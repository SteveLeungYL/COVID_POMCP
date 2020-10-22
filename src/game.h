#ifndef GAME_H
#define GAME_H

#include "simulator.h"
#include<vector>
#include<map>

using namespace std;

class GAME_STATE : public STATE
{
public:
    ///homogenous attacker
    double lambda;
};

class GAME : public SIMULATOR
{
public:
    ///constructor
    GAME(int nTarget, int nResource, double *advPenalty, double *advReward, double *defPenalty, double *defReward, map<int, vector<double> * > *MSM, map<vector<double> *, int> *inverseMSM, double MSDiscretization, double lambdaDiscretization, double maxLambdaValue);

    virtual STATE* CreateStartState() const;

    // Free memory for state
    virtual void FreeState(STATE* state) const;

    // Update state according to action, and get observation and reward. 
    // Return value of true indicates termination of episode (if episodic)
    virtual bool Step(STATE& state, int action,
        int& observation, double& reward) const;

    // Create new state and copy argument (must be same type)
    virtual STATE* Copy(const STATE& state) const;

    // Sanity check
    virtual void Validate(const STATE& state) const;

    // Modify state stochastically to some related state
    virtual bool LocalMove(STATE& state, const HISTORY& history,
        int stepObs, const STATUS& status) const;

  
    int SelectRandom(const STATE& state, const HISTORY& history,
        const STATUS& status) const;

    // Generate set of legal actions
    virtual void GenerateLegal(const STATE& state, const HISTORY& history,
        std::vector<int>& actions, const STATUS& status) const;

	 virtual void DisplayBeliefs(const BELIEF_STATE& beliefState,
        std::ostream& ostr) const;
    virtual void DisplayState(const STATE& state, std::ostream& ostr) const;
    virtual void DisplayAction(int action, std::ostream& ostr) const;
    virtual void DisplayObservation(const STATE& state, int observation, std::ostream& ostr) const;


   private:
   mutable MEMORY_POOL<GAME_STATE> MemoryPool;
   ////user defined variables
	double K;///the denominator of the MS space
   double *lambdaSpace;
   int numLambda;///the length of lambdaSpace
   double maxLambda;
   double discretization;//the discretization level of lambda
   int nTargets,nResources;
   double *advPenalties;
   double *defPenalties;
   double *advRewards;
   double *defRewards;
   map<int, vector<double> * > *MSMap;
	map<vector<double> *, int> *inverseMSMap;
};

#endif

