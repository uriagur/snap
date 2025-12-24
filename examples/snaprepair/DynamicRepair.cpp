#include "stdafx.h"
#include "ncp.h"
#include "ClusterManager.hpp"
#include "RepairStrategies.h"

// Helper function to extract nodes from local cluster result
void GetNIdV(TLocClust& Clust, int idx, TIntV& NodesInS) {
    NodesInS.Clr();
    const TIntV& AllNodes = Clust.GetNIdV();
    for (int i = 0; i <= idx && i < AllNodes.Len(); i++) {
        NodesInS.Add(AllNodes[i]);
    }
}

// Helper function to get the index of the best cut (minimum conductance)
int GetBestCutIdx(TLocClust& Clust) {
    return Clust.BestCut();
}

// The "State-of-the-Art" Dynamic Repair Implementation
// This integrates the local APPR check directly into the repair loop
void RunDynamicRepair(PUNGraph& Graph, ClusterPriorityQueue& TriageQueue, 
                     int Budget, double Threshold) {
    int totalEdgesAdded = 0;
    TLocClust Clust(Graph, 0.001); // Local worker with alpha = 0.001
    
    printf("\nRunning Dynamic Repair\n");
    printf("  Budget: %d edges\n", Budget);
    printf("  Threshold: %.3f\n", Threshold);
    
    int clustersProcessed = 0;
    int clustersFixed = 0;
    int clustersRecycled = 0;

    while (Budget > 0 && !TriageQueue.IsEmpty()) {
        ClusterInfo ci = TriageQueue.Pop();
        clustersProcessed++;

        // Progress reporting
        if (clustersProcessed % 10 == 0) {
            printf("  Progress: %d edges added, %d clusters processed, Queue size: %d\n", 
                   totalEdgesAdded, clustersProcessed, TriageQueue.Size());
        }

        // STEP 1: VERIFICATION
        // Re-run local sweep ONLY for this seed to see if it's still a "bad cut"
        Clust.FindBestCut(ci.SeedNId, ci.K, 0.001);
        int bestIdx = GetBestCutIdx(Clust);
        double currentPhi = Clust.GetPhi(bestIdx);

        // STEP 2: ORTHOGONAL REPAIR
        if (currentPhi < Threshold) {
            TIntV NodesInS;
            GetNIdV(Clust, bestIdx, NodesInS);
            
            int target = SelectTarget(Graph, NodesInS, STRATEGY_RANDOM);
            
            if (target != -1 && !Graph->IsEdge(ci.SeedNId, target) && ci.SeedNId != target) {
                Graph->AddEdge(ci.SeedNId, target);
                Budget--;
                totalEdgesAdded++;

                // STEP 3: RE-EVALUATE AND RECYCLE
                // Immediately check the new conductance after adding the edge
                Clust.FindBestCut(ci.SeedNId, ci.K, 0.001);
                int newBestIdx = GetBestCutIdx(Clust);
                double newPhi = Clust.GetPhi(newBestIdx);
                
                if (newPhi < Threshold) {
                    // Still bad, keep fixing
                    ci.Conductance = newPhi;
                    ci.LastUpdated = totalEdgesAdded;
                    TriageQueue.Push(ci);
                    clustersRecycled++;
                } else {
                    // Fixed!
                    clustersFixed++;
                }
            }
        } else {
            // Cluster is already fixed (perhaps as a side effect of fixing another cluster)
            clustersFixed++;
        }
        // If currentPhi >= Threshold, the cluster is fixed! We don't push it back.
    }
    
    printf("  Total edges added: %d\n", totalEdgesAdded);
    printf("  Clusters processed: %d\n", clustersProcessed);
    printf("  Clusters fixed: %d\n", clustersFixed);
    printf("  Clusters recycled: %d\n", clustersRecycled);
    printf("  Final queue size: %d\n", TriageQueue.Size());
}
