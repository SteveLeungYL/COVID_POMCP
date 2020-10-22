#include "game.h"
#include<time.h>
#include<math.h>
#include<string>
#include<iostream>
#include<fstream>

using namespace std;

GAME::GAME(int nTarget, int nResource, double *advPenalty, double *advReward, double *defPenalty, double *defReward, map<int, vector<double> * > *MSM, map<vector<double> *, int> *inverseMSM, double MSDiscretization, double lambdaDiscretization, double maxLambdaValue)
{
	srand(time(NULL));

	nTargets = nTarget;
	nResources = nResource;
	advPenalties = advPenalty;
	advRewards = advReward;
	defPenalties = defPenalty;
	defRewards = defReward;
	MSMap = MSM;
	inverseMSMap = inverseMSM;
	
	K = MSDiscretization;
	discretization = lambdaDiscretization;
	maxLambda = maxLambdaValue;

	lambdaSpace = new double[(int)(maxLambda/discretization)+1];
	numLambda = (int)(maxLambda/discretization)+1;//the length of lambdaSpace
	int count=0;
	for (double i=0;i<=maxLambda;i=i+discretization,count++)
		lambdaSpace[count] = i;


	/////variables from superclass
	NumActions = MSMap->size();
	NumObservations = nTargets;
	Discount = 1.0;
	RewardRange = 0.1;////redundant variable as long as autoexploration is false
}

STATE* GAME::CreateStartState() const
{
	GAME_STATE* state = MemoryPool.Allocate();
	int randomIndex = rand()%numLambda;
	state->lambda = lambdaSpace[randomIndex];
	return state;
}

void GAME::FreeState(STATE *state) const
{
   GAME_STATE* bsstate = safe_cast<GAME_STATE*>(state);
   MemoryPool.Free(bsstate);
}

STATE* GAME::Copy(const STATE& state) const
{
    assert(state.IsAllocated());
    const GAME_STATE& oldstate = safe_cast<const GAME_STATE&>(state);
    GAME_STATE* newstate = MemoryPool.Allocate();
    *newstate = oldstate;
    return newstate;
}

void GAME::Validate(const STATE& state) const
{
	bool isValid=false;
   const GAME_STATE& game_state = safe_cast<const GAME_STATE&>(state);
   for (int i=0;i<numLambda;i++)
		if (game_state.lambda==lambdaSpace[i])
		{
			isValid=true;
			break;
		}

	if (!isValid)
		assert(false);
}

void GAME::DisplayBeliefs(const BELIEF_STATE& beliefState,
    ostream& ostr) const
{
	const GAME_STATE *bsstate;
   
	ostr<<"Belief: \n";
   for (int i=0;i<beliefState.GetNumSamples();i++)
   {
      bsstate = safe_cast<const GAME_STATE*>(beliefState.GetSample(i));
   	ostr<<"State: "<<bsstate->lambda<<"\n";
   }
}

void GAME::DisplayState(const STATE& state, ostream& ostr) const
{
	const GAME_STATE& game_state = safe_cast<const GAME_STATE&>(state);
   ostr<<"State: "<<game_state.lambda<<"\n";
}

void GAME::DisplayObservation(const STATE& state, int observation, ostream& ostr) const
{
   ostr<<"Observation: Attack on Target "<<observation<<"\n";
}

void GAME::DisplayAction(int action, ostream& ostr) const
{
   ostr<<"Action: \n";
   vector<double> *observ = MSMap->find(action)->second;
   for (int i=0;i<observ->size();i++)
      ostr<<((*observ)[i])/K<<",";
   ostr<<"\n";
}

bool GAME::LocalMove(STATE& state, const HISTORY& history, int stepObs, const STATUS& status) const
{
	///we do nothing assuming that a state is stoschatically related to another state if the 2nd state can be reached from 1st by transition
	//so we do not change state
	///however this may lead to particle deprivation
	cout<<"Readding same particle\n";	
}

void GAME::GenerateLegal(const STATE& state, const HISTORY& history, vector<int>& actions, const STATUS& status) const
{
	actions.clear();
	for (map<int, vector<double> * >::iterator it=MSMap->begin();it!=MSMap->end();it++)
	{
		actions.push_back(it->first);
	}
}

bool GAME::Step(STATE& state, int action, int& observation, double& reward) const
{
	const GAME_STATE& game_state = safe_cast<const GAME_STATE&>(state);

	
	vector<double> *MS = MSMap->find(action)->second;
	double advUtilities[nTargets];
	
	for (int i=0;i<nTargets;i++)
	{
		advUtilities[i] = (1 - ((*MS)[i])/K)*advRewards[i] + (((*MS)[i])/K)*advPenalties[i];
	}
	
	//given by QR
	////attProb contains cumulative prob
	double attProb[nTargets];
	double denom=0;
	for (int i=0;i<nTargets;i++)
	{
		denom+=exp(game_state.lambda*advUtilities[i]);
		attProb[i]=denom;
	}
		
	double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	
	for (int i=0;i<nTargets;i++)
	{
		if (r<(attProb[i]/denom))
		{
			observation=i;
			break;
		}
	}
	
	
	reward = ((*MS)[observation]/K)*defRewards[observation] + (1-((*MS)[observation]/K))*defPenalties[observation];	
}
