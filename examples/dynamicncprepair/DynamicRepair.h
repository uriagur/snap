#ifndef snap_dynamic_repair_h
#define snap_dynamic_repair_h

#include "stdafx.h"
#include "ClusterManager.hpp"

// Helper function to extract nodes from local cluster result
void GetNIdV(TLocClust& Clust, int idx, TIntV& NodesInS);

// Helper function to get the index of the best cut (minimum conductance)
int GetBestCutIdx(TLocClust& Clust);

// The "State-of-the-Art" Dynamic Repair Implementation
void RunDynamicRepair(PUNGraph& Graph, ClusterPriorityQueue& TriageQueue, 
                     int Budget, double Threshold);

#endif  // snap_dynamic_repair_h
