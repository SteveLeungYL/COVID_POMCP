#ifndef SOCIAL_H
#define SOCIAL_H

#include "simulator.h"
#include "Graph.h"
#include<set>
#include<vector>

using namespace std;

class SOCIAL_STATE : public STATE
{
public:
    vector<int> netVector;
};

class SOCIAL : public SIMULATOR
{
public:
	 ///constructor
	 SOCIAL(string graphFile, int Time, int k, int numActions, vector<int> *numer);
	 
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

	 Graph *graph;
	 ////this is how the network actually looks like..comes from createStartState
	 vector<int> *groundTruth;
	 int K;///the number of nodes picked in one go
	 int numLegalActions; ////number of legal actions to return every time getLegalAction is called
	 map<int, vector<int> > actionMap;
	 map<int, vector<int> > observationMap;

	private:
	 int comb(int n, int k) const;
	 void createActionAndObservationMap(Graph *a);

	 mutable MEMORY_POOL<SOCIAL_STATE> MemoryPool;
	 //map<int, vector<int> > actionMap;
	 set<set<int> > observationMaps;
};

#endif
