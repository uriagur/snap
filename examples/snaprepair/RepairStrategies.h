#ifndef snap_repair_strategies_h
#define snap_repair_strategies_h

#include "stdafx.h"

// Strategy types for target selection
enum StrategyType {
    STRATEGY_RANDOM,     // Boundary-biased random strategy
    STRATEGY_SPECTRAL    // Future: Spectral/Fiedler-based strategy
};

// Check if a node is in the cluster
bool IsInCluster(int nid, const TIntV& ClusterNodes);

// Identify boundary nodes in the cluster
void GetBoundaryNodes(const PUNGraph& Graph, const TIntV& CurrentClusterNodes, TIntV& BoundaryNodes);

// Select target node based on strategy
int SelectTarget(const PUNGraph& Graph, const TIntV& CurrentClusterNodes, StrategyType type);

#endif  // snap_repair_strategies_h
