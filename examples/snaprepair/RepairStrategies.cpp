#include "stdafx.h"
#include "RepairStrategies.h"

// Check if a node is in the cluster
bool IsInCluster(int nid, const TIntV& ClusterNodes) {
    for (int i = 0; i < ClusterNodes.Len(); i++) {
        if (ClusterNodes[i] == nid) {
            return true;
        }
    }
    return false;
}

// Identify boundary nodes in the cluster
// Boundary nodes are nodes in S with at least one neighbor in V \ S
void GetBoundaryNodes(const PUNGraph& Graph, const TIntV& CurrentClusterNodes, TIntV& BoundaryNodes) {
    BoundaryNodes.Clr();
    
    // Create a set for fast lookup
    TIntSet ClusterSet(CurrentClusterNodes.Len());
    for (int i = 0; i < CurrentClusterNodes.Len(); i++) {
        ClusterSet.AddKey(CurrentClusterNodes[i]);
    }
    
    // Check each node in cluster
    for (int i = 0; i < CurrentClusterNodes.Len(); i++) {
        int nid = CurrentClusterNodes[i];
        if (!Graph->IsNode(nid)) continue;
        
        TUNGraph::TNodeI NI = Graph->GetNI(nid);
        // Check if node has any neighbor outside the cluster
        for (int e = 0; e < NI.GetDeg(); e++) {
            int nbr = NI.GetNbrNId(e);
            if (!ClusterSet.IsKey(nbr)) {
                BoundaryNodes.Add(nid);
                break;  // Found external neighbor, this is a boundary node
            }
        }
    }
}

// Select target node based on strategy
// STRATEGY_RANDOM: Boundary-biased random strategy
int SelectTarget(const PUNGraph& Graph, const TIntV& CurrentClusterNodes, StrategyType type) {
    if (type == STRATEGY_RANDOM) {
        // 1. Identify the boundary nodes in the cluster
        TIntV BoundaryNodes;
        GetBoundaryNodes(Graph, CurrentClusterNodes, BoundaryNodes);
        
        // If no boundary nodes, use any node from cluster
        if (BoundaryNodes.Len() == 0) {
            if (CurrentClusterNodes.Len() == 0) {
                return -1;
            }
            BoundaryNodes = CurrentClusterNodes;
        }
        
        // 2. Pick a random node from the boundary to maximize expansion
        int u = BoundaryNodes[TInt::Rnd.GetUniDevInt(BoundaryNodes.Len())];
        
        // 3. Pick a random node from the entire graph (Far-field approximation)
        if (Graph->GetNodes() == 0) {
            return -1;
        }
        
        int v = Graph->GetRndNId();
        int maxAttempts = 100;
        int attempts = 0;
        
        // Try to find a node not in cluster
        while (IsInCluster(v, CurrentClusterNodes) && attempts < maxAttempts) {
            v = Graph->GetRndNId();
            attempts++;
        }
        
        // If we couldn't find a node outside cluster after many attempts, 
        // just return v anyway (edge case for very small graphs)
        return v;
    }
    
    // Future: Add Strategy A (Spectral/Fiedler-based) here
    // if (type == STRATEGY_SPECTRAL) { ... }
    
    return -1;
}
