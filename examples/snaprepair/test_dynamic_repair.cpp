#include "stdafx.h"
#include "ncp.h"
#include "ClusterManager.hpp"
#include "RepairStrategies.h"
#include "DynamicRepair.h"

// Build Triage Queue for Dynamic Repair
void BuildTriageQueueForDynamicRepair(const PUNGraph& Graph, 
                                      ClusterPriorityQueue& TriageQueue,
                                      int KMin, int KMax, double KFac, 
                                      int Coverage, double PhiThreshold) {
    printf("\nBuilding Triage Queue for Dynamic Repair\n");
    printf("  KMin=%d, KMax=%d, KFac=%.2f, Coverage=%d, PhiThreshold=%.3f\n", 
           KMin, KMax, KFac, Coverage, PhiThreshold);
    
    TLocClust Clust(Graph, 0.001);  // Alpha = 0.001
    int totalScans = 0;
    int addedToQueue = 0;
    
    int prevK = -1;
    for (int K = KMin; K < KMax; K = int(KFac * double(K)) + 1) {
        if (K == prevK) { K++; }
        prevK = K;
        
        const double Eps = 1.0 / (10.0 * K);
        const int Runs = 2 + int(Coverage * floor(double(Graph->GetEdges()) / double(K)));
        
        printf("  K=%d, Eps=%.6f, Runs=%d\n", K, Eps, Runs);
        
        for (int run = 0; run < Runs; run++) {
            const int SeedNId = Graph->GetRndNId();
            totalScans++;
            
            // Run FindBestCut which does ApproxPageRank, normalization, and sweep
            Clust.FindBestCut(SeedNId, K, 0.001);
            
            int bestIdx = Clust.BestCut();
            double phi = Clust.GetPhi(bestIdx);
            
            if (phi < PhiThreshold) {
                TriageQueue.Push(ClusterInfo(phi, SeedNId, K, 0));
                addedToQueue++;
            }
        }
    }
    
    printf("  Total scans: %d\n", totalScans);
    printf("  Clusters added to queue: %d\n", addedToQueue);
    printf("  Queue size: %d\n", TriageQueue.Size());
}

int main(int argc, char* argv[]) {
    Env = TEnv(argc, argv, TNotify::StdNotify);
    Env.PrepArgs(TStr::Fmt("Dynamic Repair Test. build: %s, %s. Time: %s", 
                           __TIME__, __DATE__, TExeTm::GetCurTm()));
    TExeTm ExeTm;
    Try
    
    // Command line arguments
    const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", 
                                             "Input undirected graph (one edge per line)");
    const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "dynamic_repair", 
                                              "Output file prefix");
    const int Budget = Env.GetIfArgPrefixInt("-b:", 100, 
                                             "Edge budget (number of edges to add)");
    const double PhiThreshold = Env.GetIfArgPrefixFlt("-t:", 0.2, 
                                                       "Conductance threshold for valley detection");
    const int KMin = Env.GetIfArgPrefixInt("-kmin:", 10, 
                                          "Minimum cluster size K");
    const int KMax = Env.GetIfArgPrefixInt("-kmax:", 100, 
                                          "Maximum cluster size K");
    const double KFac = Env.GetIfArgPrefixFlt("-kfac:", 1.2, 
                                             "K growth factor");
    const int Coverage = Env.GetIfArgPrefixInt("-c:", 10, 
                                              "Coverage (scans per node)");
    
    printf("Dynamic Repair Test\n");
    printf("===================\n");
    printf("Input graph: %s\n", InFNm.CStr());
    printf("Output prefix: %s\n", OutFNm.CStr());
    printf("Budget: %d edges\n", Budget);
    printf("Phi threshold: %.3f\n", PhiThreshold);
    printf("\n");
    
    // Load graph
    printf("Loading graph...\n");
    PUNGraph Graph = TSnap::GetMxWcc(TSnap::LoadEdgeList<PUNGraph>(InFNm, 0, 1));
    printf("  Nodes: %d, Edges: %d\n", Graph->GetNodes(), Graph->GetEdges());
    
    // Build Triage Queue
    ClusterPriorityQueue TriageQueue;
    BuildTriageQueueForDynamicRepair(Graph, TriageQueue, KMin, KMax, KFac, Coverage, PhiThreshold);
    
    // Run Dynamic Repair
    RunDynamicRepair(Graph, TriageQueue, Budget, PhiThreshold);
    
    // Save results
    printf("\nSaving results...\n");
    
    // Save modified graph
    TStr GraphFNm = OutFNm + "_graph.txt";
    TSnap::SaveEdgeList(Graph, GraphFNm.CStr());
    printf("  Modified graph saved to: %s\n", GraphFNm.CStr());
    printf("  Final graph: Nodes=%d, Edges=%d\n", Graph->GetNodes(), Graph->GetEdges());
    
    Catch
    printf("\nTotal run time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
    
    return 0;
}
