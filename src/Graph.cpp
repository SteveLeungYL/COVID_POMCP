//
//  Graph.cpp
//  POMDP
//
//  Created by Amulya Yadav on 3/13/14.
//  Copyright (c) 2014 Amulya Yadav. All rights reserved.
//

#include "Graph.h"
#include<fstream>

using namespace std;

int Graph::getIndexOfNode(int node)
{
    int ind=0;
    for (map<int, map<int, pair<int, int> > >::iterator it= adjList.begin();it!=adjList.end();it++, ind++)
        if (it->first==node)
            return ind;
    
    return -1;
}


Graph::Graph(int nodes, int edges)
{
        numNodes = nodes;
        numEdges = edges;
}

Graph::Graph()
{
    
}
    
Graph::Graph(string graphFile, int time)
{
        T= time;
        int uncertainEdgeListCounter = 0;
        numUncertainEdges = 0;
        ifstream file(graphFile.c_str());
        int startVertex, endVertex, edgeCertain, edgeUncertain;
        if (file.is_open())
        {
            while(file>>startVertex>>endVertex>>edgeCertain>>edgeUncertain)
            {
                if (adjList.find(startVertex)==adjList.end())
                {
                    map <int, pair<int, int> > *vertexList = new map<int, pair<int, int> >();
                    adjList.insert(pair<int, map<int, pair<int, int> > >(startVertex, *vertexList));
                }
                
                if (adjList.find(endVertex)==adjList.end())
                {
                    map <int, pair<int, int> >*vertexList = new map<int, pair<int, int> >();
                    adjList.insert(pair<int, map<int, pair<int, int> > >(endVertex, *vertexList));
                }
                
                map<int, pair<int, int> > *list = &adjList.find(startVertex)->second;
                if (edgeCertain==0)///uncertain edge
                {
                    list->insert(pair<int, pair<int, int> >(endVertex, pair<int, int>(uncertainEdgeListCounter, edgeCertain)));
                    uncertainEdgeList.insert(pair<int, pair<int, int> >(uncertainEdgeListCounter++, pair<int, int>(startVertex, endVertex)));
                }
                else ///certain edge...will have a -1
                    list->insert(pair<int, pair<int, int> >(endVertex, pair<int, int>(-1, edgeCertain)));
                
                numUncertainEdges+=edgeUncertain;
                numEdges+=(edgeCertain+edgeUncertain);
            }
            file.close();
            numNodes = (int)adjList.size();
            
            ////initialize adjMatrix
            adjMatrix = new double*[numNodes];
            
            for (int i=0;i<numNodes;i++)
            {
                adjMatrix[i] = new double[numNodes];
            }
            
            int nodeIndex=0;
            for (map<int, map<int, pair<int, int> > >::iterator it=adjList.begin();it!=adjList.end();it++, nodeIndex++)
            {
                map<int, pair<int, int> > curra = it->second;
                for (map<int, pair<int, int> >::iterator it1=curra.begin();it1!=curra.end();it1++)
                {
                    if (it1->second.second==1)
                        adjMatrix[nodeIndex][getIndexOfNode(it1->first)]=1;
                    else///if uncertain
                        adjMatrix[nodeIndex][getIndexOfNode(it1->first)]=0.5;
                }
                
                ////fill up rest of entries in curr row with 0
                for (int i=0;i<numNodes;i++)
                {
                    if (adjMatrix[nodeIndex][i]!=1&&adjMatrix[nodeIndex][i]!=0.5)
                        adjMatrix[nodeIndex][i]=0;
                }
            }
            ////adjMatrix initialized
            
            ////calculate diffusion centrality vector
            
            double temp[numNodes][numNodes];
            double temp2[numNodes][numNodes];
            
            diffusionCentrality = new double[numNodes];

				columnDiffusionCentrality = new double[numNodes];
            
				///this variable is just to make sure that both colDiff and diffusionCentr can be populated in the same set of loops
				int altIndex=0;
            for (int i=0;i<numNodes;i++)
            {
                diffusionCentrality[i]=0;
					 columnDiffusionCentrality[i]=0;
                for (int j=0;j<numNodes;j++)
                {
                    temp[i][j]=adjMatrix[i][j];
                    temp2[i][j]=adjMatrix[i][j];
                    diffusionCentrality[i] += 0.5*temp[i][j];
						  columnDiffusionCentrality[i] += 0.5*adjMatrix[j][altIndex];
                }
					 altIndex++;
            }
            
            
            double prob = 0.25;
            for (int i=2;i<=T;i++, prob = prob*0.5)
            {
                for (int j=0;j<numNodes;j++)
                {
                    for (int k=0;k<numNodes;k++)
                    {
								double sum=0;
                        for (int h=0;h<numNodes;h++)
                            sum+=temp2[j][h]*adjMatrix[h][k];
                        
                        temp[j][k]=sum;
                        diffusionCentrality[j] += prob*temp[j][k];
                    }
                }
	
					for (int j=0;j<numNodes;j++)
					{
						for (int k=0;k<numNodes;k++)
						{
							columnDiffusionCentrality[j] += prob*temp[k][j];
                	}
					}


                for (int j=0;j<numNodes;j++)
                    for (int k=0;k<numNodes;k++)
                        temp2[j][k]=temp[j][k];
            }
            
            ////diffusion centrality vector made
            
            cout<<"DIFFUSION CENTRALITY\n";
            for (int h=0;h<numNodes;h++)
                cout<<diffusionCentrality[h]<<" ";
            
        }
        else
        {
            cout<<"Error in opening file in Graph constructor function";
        }
}
