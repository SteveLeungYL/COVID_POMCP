//
//  Graph.h
//  POMDP
//
//  Created by Amulya Yadav on 3/13/14.
//  Copyright (c) 2014 Amulya Yadav. All rights reserved.
//

#ifndef __POMDP__Graph__
#define __POMDP__Graph__

#include <iostream>
#include<map>
#include<string>

using namespace std;

class Graph{
public:
    int numNodes;
    int numEdges;
    int numUncertainEdges;
    
    int T; ////the time t used in the diffusion centrality calculation
    
    double **adjMatrix;
    
    double *diffusionCentrality;
	 
	 ////diffusionCentrality[i] measures the influence spreading power of ith node.
	 ////in order to calculate diffusionCentrality[i], we sum up the rows.
    ////if instead we sum up the columns, we would get a measure of the the power
	 /////of every node getting influenced...as we would be counting the number of 
	 ////paths of length k(k=1..T) that end at the ith node.
    /////thus...columnDiffusionCentrality would have the power of every node to get influenced.
    ////when you want a black box simulator for the POMDP, and you have to come up with a new state
    ///you can sample according to the probabiltities(normalized) in columnDiffusionCentrality to decide
	 ////whether to include or remove every node 
	 double *columnDiffusionCentrality;
    
    /////increment edge number only in case there is an uncertain edge
    map<int, map<int, pair<int, int> > > adjList;
    
    map<int, pair<int, int> > uncertainEdgeList;
    
    Graph(int nodes, int edges);
    
    Graph(string graphFile, int t);
    
    Graph();
    
    int getIndexOfNode(int node);
};

#endif /* defined(__POMDP__Graph__) */
