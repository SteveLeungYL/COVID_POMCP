#include "social.h"
#include<time.h>
#include<math.h>
using namespace std;

int SOCIAL::comb(int n, int k) const
{
	if (n<k)
		return 0;
	if (n==k)
		return 1;
	
    int i,j;
    int a[n][n];
    for(i=0;i<=n;i++)
    {
       for(j=0;j<=(i>k?k:i);j++)
       {
            if(j==0||j==i)  //if j==0 or j==i , coefficient =1 as nC0=1 and nCn=1 where C= combination
            {
                a[i][j]=1;
            }
            else
                a[i][j]=a[i-1][j-1]+a[i-1][j];  // nCk=(n-1)C(k-1) + (n-1)C(k)
        }
	 }
    return a[n][k];// return value of nCk
}

////INCLUDES BOTH ACTION AND OBSERVATION PROB
void SOCIAL::createActionAndObservationMap(Graph *a)
{
    ///enter null observation in map
    //vector<int> *nullObs = new vector<int>();
    
    //observationMap.insert(pair<unsigned long, vector<int> >(0, *nullObs));
	 //map<int, vector<int> > actionMap;
    
    
    int actionIndex=0;
    int observationIndex = 0;
    vector<int> tempvec;
    
    for (int i=0;i<a->numNodes;i++)
    {
        if (i>a->numNodes-1-K)
            tempvec.push_back(1);
        else
            tempvec.push_back(0);
    }
    
    bool end=false;
    
    while(!end)
    {
        vector<int> *newvec = new vector<int>();
        for (int i=0;i<tempvec.size();i++)
            newvec->push_back(tempvec[i]);
        actionMap.insert(pair<int, vector<int> >(actionIndex++, *newvec));
        
        /////get the observation associated with this action
        set<int> *currentObservation = new set<int>();
        int index=0;

        for (map<int, map<int, pair<int, int> > >::iterator it=a->adjList.begin();it!=a->adjList.end();it++, index++)
        {
            if (tempvec[index]==1)
            {                map<int, pair<int, int> > currList = it->second;
                for (map<int, pair<int, int> >::iterator listIt = currList.begin(); listIt!=currList.end();listIt++)
                {
                    if ((listIt->second).second==0)///uncertain edge
                        currentObservation->insert((listIt->second).first);///these are the indexes of the uncertainEdgeList
                    else
                        continue;
                }
            }
        }
        
        if (observationMaps.find(*currentObservation)==observationMaps.end())
        {
            ////this particular observation set has not been seen before...add all permutations of it in observationMap
            observationMaps.insert(*currentObservation);
            
            ////generate permutations of this observation
            ///-1 denotes that that particular edge was not observed during the observation..0 denotes that edge not present..1 present
            
            vector<int> tempvector;
            for (int i=0;i<currentObservation->size();i++)
                tempvector.push_back(0);
            
            bool finish=false;
            
            while(!finish)
            {
                vector<int> *newObservationVector = new vector<int>();
                for (int i=0;i<a->uncertainEdgeList.size();i++)
                    newObservationVector->push_back(-1);
                
                int permutIndex = 0;
                for (set<int>::iterator it = currentObservation->begin();it!=currentObservation->end();it++, permutIndex++)
                {
                    newObservationVector->at(*it) = tempvector[permutIndex];
                }
                
                
                observationMap.insert(pair<int, vector<int> >(observationIndex++, *newObservationVector));
                
                ////add 1 to tempvec to get new tempvec
                int index = tempvector.size()-1;
                while(index>=0 && tempvector[index]==1)
                    tempvector[index--]=0;
                
                if (index<0)
                {
                    finish=true;
                    continue;
                }
                else
                    tempvector[index]=1;
                
            }
            
            
        }/////observations added for current Action..if this if condition is bypassed that means that observation was already accounted for
        
        ////generate next binary string with k bits set
        bool nextPresent=false;
        for (int i=tempvec.size()-1;i>=1;i--)
        {
            if (tempvec[i]==1&&tempvec[i-1]==0)
            {
                nextPresent=true;
                tempvec[i]=0;
                tempvec[i-1]=1;
                int numOnesAfterPoint = 0;
                for (int j=i+1;j<tempvec.size();j++)
                    if (tempvec[j]==1)
                        numOnesAfterPoint++;
                
                for (int j=tempvec.size()-1;j>=i+1;j--)
                {
                    if (numOnesAfterPoint)
                    {
                        tempvec[j]=1;
                        numOnesAfterPoint--;
                    }
                    else
                        tempvec[j]=0;
                }
                break;
            }
            else
                continue;
        }
        
        if (!nextPresent)
            end=true;
        
    }
	
	////freeing up the action map
	/*cout<<"Action map: "<<actionMap.size()<<"\n";
	for (map<int, vector<int> >::iterator ww=actionMap.begin();ww!=actionMap.end();ww++)
	{
		vector<int> d=ww->second;
		for (int y=0;y<d.size();y++)
			cout<<d[y]<<",";
		cout<<"\n";
	//	delete &(ww->second);
	}*/
}


SOCIAL::SOCIAL(string graphFile, int Time, int k, int numActions, vector<int>* rec)
{
	graph = new Graph(graphFile, Time);
	K=k;
	groundTruth = new vector<int>();
	
	for (int i=0;i<rec->size();i++)
      groundTruth->push_back(rec->at(i));
	
   numLegalActions = numActions;
   //NumObservations = NumActions; ////an upper bound see notes for explanation
   Discount = 0.95;
   RewardRange = 0.1;
   srand(time(NULL));
	createActionAndObservationMap(graph);
	NumActions = actionMap.size();
   NumObservations = observationMap.size();
	//cout<<"Number of actions in constructor: "<<actionMap.size()<<"\n";
   //cout<<"Number of observations in constructor: "<<observationMap.size()<<"\n";
	//cout<<"Number of actions in constructor: "<<NumActions<<"\n";
   //cout<<"Number of observations in constructor: "<<NumObservations<<"\n";
	cout<<"Number of actions in constructor: "<<GetNumActions()<<"\n";
	cout<<"Number of observations in constructor: "<<GetNumObservations()<<"\n";
}

STATE* SOCIAL::CreateStartState() const
{
	SOCIAL_STATE* bsstate = MemoryPool.Allocate();

	/////initially all nodes are uninfluenced
	for (int i=0;i<graph->numNodes;i++)
		bsstate->netVector.push_back(0);

   ////randomly set every uncertain edge to either 0 or 1
	///this is the way it has been done in battleship
	///the only care taken is to make the state valid
	for (int i=0;i<graph->numUncertainEdges;i++)
		bsstate->netVector.push_back(rand()%2);
	
//	cout<<"netvector size: "<<bsstate->netVector.size()<<"\n";
//	cout<<"groundTruth size: "<<groundTruth->size()<<"\n";

	///copy startstate into groundTruth
	/*for (int i=0;i<bsstate->netVector.size();i++)
	{
		//cout<<"I val: "<<i<<"\n";
		groundTruth->push_back(bsstate->netVector[i]);
	}*/
	//cout<<"Value returned\n";
	return bsstate;
}

void SOCIAL::FreeState(STATE *state) const
{
	SOCIAL_STATE* bsstate = safe_cast<SOCIAL_STATE*>(state);
   MemoryPool.Free(bsstate);
}

STATE* SOCIAL::Copy(const STATE& state) const
{
	 assert(state.IsAllocated());
    const SOCIAL_STATE& oldstate = safe_cast<const SOCIAL_STATE&>(state);
    SOCIAL_STATE* newstate = MemoryPool.Allocate();
    *newstate = oldstate;
    return newstate;
}

void SOCIAL::Validate(const STATE& state) const
{
	const SOCIAL_STATE& bsstate = safe_cast<const SOCIAL_STATE&>(state);
	for (int i=0;i<bsstate.netVector.size();i++)
	{
		if (bsstate.netVector[i]==1||bsstate.netVector[i]==0)
			continue;
		else
			assert(false);
	}
}

void SOCIAL::DisplayBeliefs(const BELIEF_STATE& beliefState,
    ostream& ostr) const
{
	cout<<"Belief: \n";
	for (int i=0;i<beliefState.GetNumSamples();i++)
	{
		cout<<"State: \n";
		const SOCIAL_STATE *bsstate = safe_cast<const SOCIAL_STATE*>(beliefState.GetSample(i));
		//const SOCIAL_STATE& bsstate = safe_cast<SOCIAL_STATE&>(beliefState.GetSample(i));
		 for (int i=0;i<bsstate->netVector.size();i++)
      	cout<<bsstate->netVector[i]<<",";
   	 cout<<"\n";
	}
}

void SOCIAL::DisplayState(const STATE& state, ostream& ostr) const
{
	cout<<"State: \n";
	const SOCIAL_STATE& bsstate = safe_cast<const SOCIAL_STATE&>(state);
	for (int i=0;i<bsstate.netVector.size();i++)
		cout<<bsstate.netVector[i]<<",";
	cout<<"\n";
}

void SOCIAL::DisplayObservation(const STATE& state, int observation, ostream& ostr) const
{
	cout<<"Observation: "<<observation<<"\n";
	/*vector<int> num;
	while(observation>0)
	{
		num.insert(num.begin(),observation%3);
		observation=observation/3;
	}

	 while(num.size()<graph->numUncertainEdges)
		num.insert(num.begin(),0);
	
	for (int i=0;i<num.size();i++)
		cout<<(num[i]-1)<<",";
	cout<<"\n";*/

	vector<int> observ = observationMap.find(observation)->second;
	for (int i=0;i<observ.size();i++)
		cout<<observ[i]<<",";
	cout<<"\n";
}

void SOCIAL::DisplayAction(int action, ostream& ostr) const
{
	cout<<"Action: \n";
	/*vector<int> index;
   ////find the corresponding action
   int elemFound=0, choose;
   int RHS = action;
   int val=K;

   while(elemFound<K)
   {
      choose = val-1;
      ///find elemFound element
      while(comb(choose,val)<RHS)
      {
         choose++;
      }
      index.push_back(choose-1);
      RHS = RHS - comb(choose-1, val);
      val=val-1;
      elemFound++;
   }

   vector<int> actionVec;
   int currIndex=0;
   for (int i=index.size()-1;i>=0;i--)
   {
      int num=index[i];
      for (int j=0;j<num-currIndex;j++)
         actionVec.push_back(0);
      actionVec.push_back(1);
      currIndex=num+1;
   }

   while(actionVec.size()<graph->numNodes)
      actionVec.push_back(0);

	for (int i=0;i<actionVec.size();i++)
		cout<<actionVec[i]<<",";
	cout<<"\n";*/

	vector<int> observ = actionMap.find(action)->second;
   for (int i=0;i<observ.size();i++)
      cout<<observ[i]<<",";
   cout<<"\n";
}

bool SOCIAL::LocalMove(STATE& state, const HISTORY& history, int stepObs, const STATUS& status) const
{
	SOCIAL_STATE& bsstate = safe_cast<SOCIAL_STATE&>(state);
	vector<int> unionState;
	for (int i=0;i<bsstate.netVector.size();i++)
		unionState.push_back(0);
	
	for (int t=0;t<history.Size();t++)
	{
		 int action = history[t].Action;
		 vector<int> actionVec = actionMap.find(action)->second;
		 ///create actionVec;
		 /*vector<int> index;
       ////find the corresponding action
       int elemFound=0, choose;
       int RHS = action;
       int val=K;

       while(elemFound<K)
       {
       	choose = val-1;
      	///find elemFound element
     		while(comb(choose,val)<RHS)
      	{
        		choose++;
      	}
      	index.push_back(choose-1);
      	RHS = RHS - comb(choose-1, val);
      	val=val-1;
      	elemFound++;
   	} 

   	if (index.size()!=K)
      	cout<<"ERROR IN GENERATING ACTION\n";

   	vector<int> actionVec;
   	int currIndex=0;
   	for (int i=index.size()-1;i>=0;i--)
   	{
      	int num=index[i]; 
      	for (int j=0;j<num-currIndex;j++)
         	actionVec.push_back(0);
      	actionVec.push_back(1);
      	currIndex=num+1;
   	}

   	while(actionVec.size()<graph->numNodes)
      	actionVec.push_back(0);*/
		
		for (int j=0;j<graph->numNodes;j++)
		{
			if (actionVec[j]==1)
				unionState[j]=1;
		}

		int observ = history[t].Observation;

		vector<int> num = observationMap.find(observ)->second;

		/*vector<int> num;
   	while(observ>0)
   	{
      	num.insert(num.begin(),(observ%3)-1);
      	observ=observ/3;
   	}
		
		while(num.size()<graph->numUncertainEdges)
      num.insert(num.begin(),-1);*/
		
		for (int j=graph->numNodes;j<bsstate.netVector.size();j++)
		{
			if (num[j-graph->numNodes]==0||num[j-graph->numNodes]==1)
				unionState[j]=1;///you've seen this edge
		}
	}

	for (int i=0;i<bsstate.netVector.size();i++)
	{
		if (unionState[i]==0)
			bsstate.netVector[i]=rand()%2;
	}
		
	return true;
}

int SOCIAL::SelectRandom(const STATE& state, const HISTORY& history, const STATUS& status) const
{
	const SOCIAL_STATE& bsstate = safe_cast<const SOCIAL_STATE&>(state);
	
	vector<int> unionState;
	for (int i=0;i<graph->numNodes;i++)
		unionState.push_back(0);
	
	for (int t=0;t<history.Size();t++)
	{
		int action = history[t].Action;
		vector<int> actionVec = actionMap.find(action)->second;
       ///create actionVec;
       /*vector<int> index;
       ////find the corresponding action
       int elemFound=0, choose;
       int RHS = action;
       int val=K;

       while(elemFound<K)
       {
         choose = val-1;
         ///find elemFound element
         while(comb(choose,val)<RHS)
         {
            choose++;
         }
         index.push_back(choose-1);
         RHS = RHS - comb(choose-1, val);
         val=val-1;
         elemFound++;
      }

      if (index.size()!=K)
         cout<<"ERROR IN GENERATING ACTION\n";

      vector<int> actionVec;
      int currIndex=0;
      for (int i=index.size()-1;i>=0;i--)
      {
         int num=index[i];
         for (int j=0;j<num-currIndex;j++)
            actionVec.push_back(0);
         actionVec.push_back(1);
         currIndex=num+1;
      }

      while(actionVec.size()<graph->numNodes)
         actionVec.push_back(0);
		*/
      for (int j=0;j<graph->numNodes;j++)
      {
         if (actionVec[j]==1)
            unionState[j]=1;
      }
	}
	
	/////generate random legal action...make sure that action picks nodes which have not been picked before
	vector<int> zeroIndices;
	vector<int> oneIndices;
	int numZeros=0;
	for (int i=0;i<unionState.size();i++)
	{
		if (unionState[i]==0)
		{
			numZeros++;
			zeroIndices.push_back(i);
		}
		else
		{
			oneIndices.push_back(i);
		}
	}

	int numChosen=0;
	
	vector<int> chosenAction;
	for (int i=0;i<graph->numNodes;i++)
		chosenAction.push_back(0);
	
	//cout<<"zeroIndices.size(): "<<zeroIndices.size()<<"\n";
	
	while(zeroIndices.size()>0&&numChosen<K)
	{
		int rando = rand()%(zeroIndices.size());
		chosenAction[zeroIndices[rando]]=1;
		zeroIndices.erase(zeroIndices.begin()+rando);
		numChosen++;
	}
	
	///only if you have less than K unpicked nodes...now need to pick from already picked nodes
	while(numChosen<K)
	{
		//cout<<"Hello\n";
		int rando = rand()%(oneIndices.size());
      chosenAction[oneIndices[rando]]=1;
      oneIndices.erase(oneIndices.begin()+rando);
		numChosen++;
	}
	
	int resultAction=0;
	///get int representation of actionVec
	/*int resultAction=0;
	int num=1;
	for (int i=0;i<chosenAction.size();i++)
	{
		if (chosenAction[i]==1)
		{
			resultAction += comb(i,num);
			num++;
		}
	}*/
	int indexaw=0;
	for (map<int, vector<int> >::const_iterator it=actionMap.begin();it!=actionMap.end();it++,indexaw++)
	{
		vector<int> currAcc = it->second;
		bool same=true;
		for (int f=0;f<currAcc.size();f++)
		{
			if (currAcc[f]!=chosenAction[f])
			{
				same=false;
				break;
			}	
		}
		if (same)
		{
			resultAction=indexaw;
			break;
		}
	}
	
	return resultAction;
}

void SOCIAL::GenerateLegal(const STATE& state, const HISTORY& history, vector<int>& legal, const STATUS& status) const
{
	const SOCIAL_STATE& bsstate = safe_cast<const SOCIAL_STATE&>(state);
   vector<int> unionState;
   for (int i=0;i<graph->numNodes;i++)
      unionState.push_back(0);

   for (int t=0;t<history.Size();t++)
   {
      int action = history[t].Action;
		vector<int> actionVec = actionMap.find(action)->second;
	   ///create actionVec;
       /*vector<int> index;
       ////find the corresponding action
       int elemFound=0, choose;
       int RHS = action;
       int val=K;

       while(elemFound<K)
       {
         choose = val-1;
         ///find elemFound element
         while(comb(choose,val)<RHS)
         {
            choose++;
         }
         index.push_back(choose-1);
         RHS = RHS - comb(choose-1, val);
         val=val-1;
         elemFound++;
      }

      if (index.size()!=K)
         cout<<"ERROR IN GENERATING ACTION\n";

      vector<int> actionVec;
      int currIndex=0;
      for (int i=index.size()-1;i>=0;i--)
      {
         int num=index[i];
         for (int j=0;j<num-currIndex;j++)
            actionVec.push_back(0);
         actionVec.push_back(1);
         currIndex=num+1;
      }

      while(actionVec.size()<graph->numNodes)
         actionVec.push_back(0);
		*/

      for (int j=0;j<graph->numNodes;j++)
      {
         if (actionVec[j]==1)
            unionState[j]=1;
      }
   }

	legal.clear();
	vector<int> zeroIndices;
   vector<int> oneIndices;

	
	int m=0;
	while(m<numLegalActions)///loop start
	{
	zeroIndices.clear();
	oneIndices.clear();
   int numZeros=0;
   for (int i=0;i<unionState.size();i++)
   {
      if (unionState[i]==0)
      {
         numZeros++;
         zeroIndices.push_back(i);
      }
      else
      {
         oneIndices.push_back(i);
      }
   }

   int numChosen=0;

   vector<int> chosenAction;
   for (int i=0;i<graph->numNodes;i++)
      chosenAction.push_back(0);

   while(zeroIndices.size()>0&&numChosen<K)
   {
      int rando = rand()%(zeroIndices.size());
      chosenAction[zeroIndices[rando]]=1;
      zeroIndices.erase(zeroIndices.begin()+rando);
      numChosen++;
   }
   
   ///only if you have less than K unpicked nodes...now need to pick from already picked nodes
   while(numChosen<K)
   {
      int rando = rand()%(oneIndices.size());
      chosenAction[oneIndices[rando]]=1;
      oneIndices.erase(oneIndices.begin()+rando);
      numChosen++;
   }

	int resultAction=0;

   ////get int representation of actionVec
   /*int resultAction=0;
   int num=1;
   for (int i=0;i<chosenAction.size();i++)
   {
		if (chosenAction[i]==1)
      {
         resultAction += comb(i,num);
         num++;
      }
   }*/

	int indexaw=0;
   for (map<int, vector<int> >::const_iterator it=actionMap.begin();it!=actionMap.end();it++,indexaw++)
   {
      vector<int> currAcc = it->second;
      bool same=true;
      for (int f=0;f<currAcc.size();f++)
      {
         if (currAcc[f]!=chosenAction[f])
         {
            same=false;
            break;
         }  
      }
      if (same)
      {
         resultAction=indexaw;
         break;
      }
   }

	///check if resultAction is already in legal...if it is...don't put it in..else put it in
	/*bool duplicate=false;
	for (int i=0;i<legal.size();i++)
	{
		if (legal[i]==resultAction)
		{
			duplicate=true;
			break;
		}
	}

	if (!duplicate)
	{*/
		legal.push_back(resultAction);
		m++;
	//}
	//else
	//	continue;

	}///loop end
	
		
}

bool SOCIAL::Step(STATE& state, int action, int& observation, double& reward) const
{
	//cout<<"Step begins\n";
   SOCIAL_STATE& bsstate = safe_cast<SOCIAL_STATE&>(state);
	
	int initialVal=0;
	for (int i=0;i<graph->numNodes;i++)
		if (bsstate.netVector[i]==1)
			initialVal++;

   vector<int> actionVec = actionMap.find(action)->second;
	
	for (int i=0;i<graph->numNodes;i++)
   {
      if (actionVec[i]==1&&bsstate.netVector[i]==1)
      {
         reward=0;
         observation=0;
         assert(false);
      }
   }

	
	double **newadjMatrix = new double*[graph->numNodes];

   for (int i=0;i<graph->numNodes;i++)
	{
		newadjMatrix[i] = new double[graph->numNodes];
		for (int j=0;j<graph->numNodes;j++)
			newadjMatrix[i][j]=0;
	}
	//cout<<"Reached uptill here\n";

	set<int> *currentObservation = new set<int>();
   int indexa=0;
   for (map<int, map<int, pair<int, int> > >::iterator it= graph->adjList.begin();it!=graph->adjList.end();it++, indexa++)
   {
      if (actionVec[indexa]==1)
      {
			//cout<<"YAHOO\n";
          map<int, pair<int, int> > currList = it->second;
          for (map<int, pair<int, int> >::iterator listIt = currList.begin(); listIt!=currList.end();listIt++)
          {
				 int endVertex=0;
             int count=0;
             for (map<int, map<int, pair<int, int> > >::iterator itee=graph->adjList.begin();itee!=graph->adjList.end();itee++,count++)
             {
					 //cout<<itee->first<<","<<listIt->first<<"\n";
                if (itee->first==listIt->first)
                {
                     endVertex = count;
							break;
                }
             }
				  ////end vertex found

             if ((listIt->second).second==0)///uncertain edge
				 {
					 //listIt->second.second=1; ////no longer an uncertain edge
					 newadjMatrix[indexa][endVertex]=groundTruth->at(graph->numNodes+(listIt->second).first);
                currentObservation->insert((listIt->second).first);///these are the indexes of the uncertainEdgeList
				 }
             else////certain edge
             {
					newadjMatrix[indexa][endVertex] = 1;
				 }
          }
      }
		else if (bsstate.netVector[indexa]==1)
		{
			//cout<<"HANNNAD\n";
			 map<int, pair<int, int> > currList = it->second;
          for (map<int, pair<int, int> >::iterator listIt = currList.begin(); listIt!=currList.end();listIt++)
          {
             int endVertex=0;
             int count=0;
             for (map<int, map<int, pair<int, int> > >::iterator itee=graph->adjList.begin();itee!=graph->adjList.end();itee++,count++)
             {
                if (itee->first==listIt->first)
                {
                     endVertex = count;
							break;
                }
             }
              ////end vertex found
				if ((listIt->second).second==0)///uncertain edge
             {  
                
                newadjMatrix[indexa][endVertex]=0.5;
             }
             else////certain edge
             {  
               newadjMatrix[indexa][endVertex] = 1;
             }
			} 
		}
		else
			continue;
   }

	//cout<<"Step constructed adjMatrix\n";
				double *diffusionCentrality = new double[graph->numNodes];
				double *columnDiffusionCentrality = new double[graph->numNodes];
				double** temp = new double*[graph->numNodes];
            for(int i = 0; i < graph->numNodes; ++i)
               temp[i] = new double[graph->numNodes];

            double** temp2 = new double*[graph->numNodes];
            for(int i = 0; i < graph->numNodes; ++i)
               temp2[i] = new double[graph->numNodes];

            ///this variable is just to make sure that both colDiff and diffusionCentr can be populated in the same set of loops
            int altIndex=0;
            for (int i=0;i<graph->numNodes;i++)
            {
                diffusionCentrality[i]=0;
                columnDiffusionCentrality[i]=0;
                for (int j=0;j<graph->numNodes;j++)
                {
                    temp[i][j]=newadjMatrix[i][j];
                    temp2[i][j]=newadjMatrix[i][j];
                    diffusionCentrality[i] += 0.5*temp[i][j];
                    columnDiffusionCentrality[i] += 0.5*newadjMatrix[j][altIndex];
                }
                altIndex++;
            }


            double prob = 0.25;
            for (int i=2;i<=graph->T;i++, prob = prob*0.5)
            {
                for (int j=0;j<graph->numNodes;j++)
                {
                    for (int k=0;k<graph->numNodes;k++)
                    {
                        double sum=0;
                        for (int h=0;h<graph->numNodes;h++)
                            sum+=temp2[j][h]*newadjMatrix[h][k];

                        temp[j][k]=sum;
                        diffusionCentrality[j] += prob*temp[j][k];
                    }
                }

               for (int j=0;j<graph->numNodes;j++)
               {
                  for (int k=0;k<graph->numNodes;k++)
						{
                     columnDiffusionCentrality[j] += prob*temp[k][j];
                  }
               }


                for (int j=0;j<graph->numNodes;j++)
                    for (int k=0;k<graph->numNodes;k++)
                        temp2[j][k]=temp[j][k];
            }

            /////delete allocated memory
            for(int i = 0; i < graph->numNodes; ++i) {
					delete [] newadjMatrix[i];
               delete [] temp[i];
               delete [] temp2[i];
            }

		
		vector<int> Observation;
   for (int i=0;i<graph->numUncertainEdges;i++)
   {
      if (currentObservation->find(i)!=currentObservation->end())
      {
         Observation.push_back(groundTruth->at(i+graph->numNodes));
         ////update state according to observation
         bsstate.netVector[graph->numNodes+i]=groundTruth->at(i+graph->numNodes);
      }
      else
         Observation.push_back(-1);
   }


	observation=0;

   int obsNumber=0;
   for (map<int, vector<int> >::const_iterator it=observationMap.begin();it!=observationMap.end();it++,obsNumber++)
   {
      vector<int> obs = it->second;
      bool same=true;
      for(int c=0;c<obs.size();c++)
      {
         if (obs[c]!=Observation[c])
         {
            same=false;
            break;
         }
      }
      if (same)
      {
         observation = obsNumber;
         break;
      }
   }

	 double total=0;
   map<int, double> sampleProb;
   for (int i=0;i<graph->numNodes;i++)
   {
      if (bsstate.netVector[i]==1||actionVec[i]==1)
         bsstate.netVector[i]=1;
      else
      {
         sampleProb.insert(pair<int, double>(i, columnDiffusionCentrality[i]));
         total+=columnDiffusionCentrality[i];
      }
   }

   for (int i=0;i<graph->numNodes;i++)
   {
      if (bsstate.netVector[i]!=1)
      {
         double prob= (sampleProb.find(i)->second)/total;
         double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
         if (r<prob)
            bsstate.netVector[i]=1;
         else
            bsstate.netVector[i]=0;
      }
   }

   /////next come up with the rewards

   reward=0;
   for (int i=0;i<graph->numNodes;i++)
      reward+=bsstate.netVector[i];

	reward = reward-initialVal;
	//cout<<"Reward: "<<reward<<"\n";

   ///check is state is terminal..return true if this is the case
   bool flag=true;
   for (int i=0;i<graph->numNodes;i++)
   {
      if (bsstate.netVector[i]!=1)
      {
         flag=false;
         break;
      }
   }
	
	delete diffusionCentrality;
	delete columnDiffusionCentrality;

	return flag;

}




	




/*bool SOCIAL::Step(STATE& state, int action, int& observation, double& reward) const
{
	SOCIAL_STATE& bsstate = safe_cast<SOCIAL_STATE&>(state);
	
	vector<int> actionVec = actionMap.find(action)->second;
	/*
	vector<int> index;
	////find the corresponding action
	int elemFound=0, choose;
	int RHS = action;
	int val=K;

	while(elemFound<K)
	{
		choose = val-1;
		///find elemFound element
		while(comb(choose,val)<RHS)
		{
			choose++;
		}
		index.push_back(choose-1);
		RHS = RHS - comb(choose-1, val);
		val=val-1;
		elemFound++;
	} 

	if (index.size()!=K)
		cout<<"ERROR IN GENERATING ACTION\n";

	vector<int> actionVec;
	int currIndex=0;
	for (int i=index.size()-1;i>=0;i--)
	{
		int num=index[i];	
		for (int j=0;j<num-currIndex;j++)
			actionVec.push_back(0);
		actionVec.push_back(1);
		currIndex=num+1;
	}

	while(actionVec.size()<graph->numNodes)
		actionVec.push_back(0);
	////action found
	/////////////////////////////////////////////////THIS NEEDS TO BE a CLOSING COMMENT....JUST SOME CHANGES IN LIEU OF NEW STEP FUNC IMPLEMENTATION
	///this part ensures that you can't pick the same nodes again
	for (int i=0;i<graph->numNodes;i++)
	{
		if (actionVec[i]==1&&bsstate.netVector[i]==1)
		{
			reward=0;
			observation=0;
			assert(false);
		}
	}

	////traverse the adjancency list to see which edges we would see
	////NOTE THAT WE ARE LOOKING AT OUTGOING EDGES FROM NODES
	set<int> *currentObservation = new set<int>();
	int indexa=0;
	for (map<int, map<int, pair<int, int> > >::iterator it= graph->adjList.begin();it!=graph->adjList.end();it++, indexa++)
	{
		if (actionVec[indexa]==1)
		{
			 map<int, pair<int, int> > currList = it->second;
          for (map<int, pair<int, int> >::iterator listIt = currList.begin(); listIt!=currList.end();listIt++)
          {
             if ((listIt->second).second==0)///uncertain edge
                currentObservation->insert((listIt->second).first);///these are the indexes of the uncertainEdgeList
             else
                continue;
          }
		}
	}

	/////note that I don't need to remove the edges in currentObservation from uncertainEdgeList(Technically,
   /////I should because all these edges are no longer uncertain as I am going to observe them now.
   /////However, even if I don't remove them, I will never be able to see the same uncertain edge again in 
   /////any observation. As I don't ever pick the same nodes again. And I will only see the same uncertain
   /////edge only if I pick the same node again. Thus, in all future observations, I will definitely see a 
   ////-1 in the edges which are right now in currentObservation. 

	vector<int> Observation;
	for (int i=0;i<graph->numUncertainEdges;i++)
	{
		if (currentObservation->find(i)!=currentObservation->end())
		{
			Observation.push_back(groundTruth->at(i+graph->numNodes));
			////update state according to observation
			bsstate.netVector[graph->numNodes+i]=groundTruth->at(i+graph->numNodes);
		}
		else
			Observation.push_back(-1);
	}

	//////need to come up with a hash map to represent this vector observation as an integer
	///add 1 to every element in observation...you now have a number expresed in ternary..convert to decimal
	/*vector<int> ternaryObservation;
	for (int i=0;i<Observation.size();i++)
		ternaryObservation.push_back(Observation[i]+1);

	observation=0;
	int base=1;
	for (int i=ternaryObservation.size()-1;i>=0;i--)
	{
		observation += ternaryObservation[i]*base;
		base=base*3;
	}
   ////////////////////////////////NNEEDS TO BE A CLOSING COMMENT IF THIS FUNCTION IS EVER USED

	observation=0;
	
	int obsNumber=0;
	for (map<int, vector<int> >::const_iterator it=observationMap.begin();it!=observationMap.end();it++,obsNumber++)
	{
		vector<int> obs = it->second;
		bool same=true;
		for(int c=0;c<obs.size();c++)
		{
			if (obs[c]!=Observation[c])
			{
				same=false;
				break;
			}
		}
		if (same)
		{
			observation = obsNumber;
			break;
		}
	}
	/////come up with the next state
	
	///first create SuA
	double total=0;
	map<int, double> sampleProb;
	for (int i=0;i<graph->numNodes;i++)
	{
		if (bsstate.netVector[i]==1||actionVec[i]==1)
			bsstate.netVector[i]=1;
		else
		{
			sampleProb.insert(pair<int, double>(i, graph->columnDiffusionCentrality[i]));
			total+=graph->columnDiffusionCentrality[i];
		}
	}
	
	for (int i=0;i<graph->numNodes;i++)
	{
		if (bsstate.netVector[i]!=1)
		{
			double prob= (sampleProb.find(i)->second)/total;
			double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
			if (r<prob)
				bsstate.netVector[i]=1;
			else
				bsstate.netVector[i]=0;
		}
	}
	
	/////next come up with the rewards

	reward=0;
	for (int i=0;i<graph->numNodes;i++)
		reward+=bsstate.netVector[i];

	///check is state is terminal..return true if this is the case
	bool flag=true;
	for (int i=0;i<graph->numNodes;i++)
	{
		if (bsstate.netVector[i]!=1)
		{
			flag=false;
			break;
		}
	}
	
	return flag;

}*/
