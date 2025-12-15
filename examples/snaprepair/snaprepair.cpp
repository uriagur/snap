#include "stdafx.h"
#include "ncp.h"

// Structure to hold cluster information for the priority queue
struct ClusterInfo {
  double Conductance;
  int SeedNId;
  int K;          // Cluster size (used to recompute cluster)
  int Volume;
  int Cut;

  ClusterInfo() : Conductance(1.0), SeedNId(-1), K(0), Volume(0), Cut(0) {}

  ClusterInfo(double phi, int seed, int k, int vol, int cut) :
    Conductance(phi), SeedNId(seed), K(k), Volume(vol), Cut(cut) {}

  // For sorting (lower conductance first)
  bool operator<(const ClusterInfo& other) const {
    return Conductance < other.Conductance;
  }
};

// Min-heap priority queue for ClusterInfo (lowest conductance has highest priority)
class ClusterPriorityQueue {
private:
  TVec<ClusterInfo> Heap;
  
  int Parent(int i) { return (i - 1) / 2; }
  int Left(int i) { return 2 * i + 1; }
  int Right(int i) { return 2 * i + 2; }
  
  void HeapifyUp(int idx) {
    while (idx > 0 && Heap[idx].Conductance < Heap[Parent(idx)].Conductance) {
      Heap.Swap(idx, Parent(idx));
      idx = Parent(idx);
    }
  }
  
  void HeapifyDown(int idx) {
    int smallest = idx;
    int left = Left(idx);
    int right = Right(idx);
    
    if (left < Heap.Len() && Heap[left].Conductance < Heap[smallest].Conductance) {
      smallest = left;
    }
    if (right < Heap.Len() && Heap[right].Conductance < Heap[smallest].Conductance) {
      smallest = right;
    }
    
    if (smallest != idx) {
      Heap.Swap(idx, smallest);
      HeapifyDown(smallest);
    }
  }
  
public:
  ClusterPriorityQueue() {}
  
  void Push(const ClusterInfo& cluster) {
    Heap.Add(cluster);
    HeapifyUp(Heap.Len() - 1);
  }
  
  ClusterInfo Pop() {
    IAssert(!IsEmpty());
    ClusterInfo top = Heap[0];
    Heap[0] = Heap[Heap.Len() - 1];
    Heap.DelLast();
    if (!IsEmpty()) {
      HeapifyDown(0);
    }
    return top;
  }
  
  bool IsEmpty() const {
    return Heap.Len() == 0;
  }
  
  int Size() const {
    return Heap.Len();
  }

  // Save priority queue to binary file (space-efficient)
  void Save(const TStr& FileName) const {
    TFOut FOut(FileName);
    FOut.Save(Heap.Len());
    for (int i = 0; i < Heap.Len(); i++) {
      const ClusterInfo& ci = Heap[i];
      FOut.Save(ci.Conductance);
      FOut.Save(ci.SeedNId);
      FOut.Save(ci.K);
      FOut.Save(ci.Volume);
      FOut.Save(ci.Cut);
    }
  }

  // Load priority queue from binary file
  void Load(const TStr& FileName) {
    TFIn FIn(FileName);
    int len;
    FIn.Load(len);
    Heap.Gen(len);
    for (int i = 0; i < len; i++) {
      ClusterInfo ci;
      FIn.Load(ci.Conductance);
      FIn.Load(ci.SeedNId);
      FIn.Load(ci.K);
      FIn.Load(ci.Volume);
      FIn.Load(ci.Cut);
      Heap.Add(ci);
    }
  }
};


// Compute Fiedler vector approximation using degree-based heuristic
// The Fiedler vector (second eigenvector of graph Laplacian) provides spectral coordinates.
// This approximation uses degree centrality normalized around the median degree.
// Nodes with above-median degree get positive values, below-median get negative values.
// This provides a reasonable spectral partition for target selection.
void ComputeFiedlerVector(const PUNGraph& Graph, TIntFltH& FiedlerVec) {
  printf("Computing Fiedler vector approximation...\n");
  FiedlerVec.Clr();
  
  const int N = Graph->GetNodes();
  if (N == 0) {
    printf("  Warning: Graph has no nodes\n");
    return;
  }
  
  // Collect all node degrees
  TIntV DegV;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    DegV.Add(NI.GetDeg());
  }
  
  // Sort to find median
  DegV.Sort();
  double median;
  if (DegV.Len() == 0) {
    median = 1.0;
  } else if (DegV.Len() % 2 == 1) {
    // Odd number of elements: take middle element
    median = double(DegV[DegV.Len() / 2].Val);
  } else {
    // Even number of elements: average of two middle elements
    int mid = DegV.Len() / 2;
    median = (double(DegV[mid - 1].Val) + double(DegV[mid].Val)) / 2.0;
  }
  
  printf("  Median degree: %.1f\n", median);
  
  // Assign Fiedler values based on deviation from median
  // This creates a spectral partition: high-degree nodes (positive) vs low-degree nodes (negative)
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    int nid = NI.GetId();
    double val = (NI.GetDeg() - median) / (median + 1.0);
    FiedlerVec.AddDat(nid, val);
  }
  
  printf("  Fiedler vector approximation computed for %d nodes\n", FiedlerVec.Len());
}

// Compute Global PageRank using SNAP's built-in implementation
void ComputeGlobalPageRank(const PUNGraph& Graph, TIntFltH& PageRankVec, double Alpha=0.85) {
  printf("Computing Global PageRank...\n");
  const double Eps = 1e-4;
  const int MaxIter = 100;
  
  // Use SNAP's built-in PageRank computation
  TSnap::GetPageRank(Graph, PageRankVec, Alpha, Eps, MaxIter);
  
  printf("  PageRank computed for %d nodes\n", PageRankVec.Len());
}

// Hash function for node set
TUInt64 HashNodeSet(const TIntV& NodeIds) {
  TIntV Sorted = NodeIds;
  Sorted.Sort();
  TUInt64 hash = 0;
  for (int i = 0; i < Sorted.Len(); i++) {
    hash = hash * 31 + TUInt64(Sorted[i].Val);
  }
  return hash;
}

// Phase II: Robust SNAP Scanner
void BuildTriageQueue(const PUNGraph& Graph, 
                     ClusterPriorityQueue& TriageQueue,
                     int KMin, int KMax, double KFac, int Coverage, double PhiThreshold) {
  printf("\nPhase II: Building Triage Queue\n");
  printf("  KMin=%d, KMax=%d, KFac=%.2f, Coverage=%d, PhiThreshold=%.3f\n", 
         KMin, KMax, KFac, Coverage, PhiThreshold);
  
  TLocClust Clust(Graph, 0.001);  // Alpha = 0.001
  THash<TUInt64, TBool> SeenClusters;
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
      
      // Valley Detection: scan for all cuts below threshold
      for (int i = 0; i < Clust.GetRndWalkSup(); i++) {
        double phi = Clust.GetPhi(i);
        
        if (phi < PhiThreshold) {
          // Extract cluster node IDs for deduplication only
          TIntV ClusterNIds;
          Clust.GetNIdV().GetSubValV(0, i, ClusterNIds);
          
          // Check if we've seen this cluster
          TUInt64 hash = HashNodeSet(ClusterNIds);
          if (!SeenClusters.IsKey(hash)) {
            SeenClusters.AddDat(hash, true);
            
            int vol = Clust.GetVol(i);
            int cut = Clust.GetCut(i);
            int clusterSize = ClusterNIds.Len();

            // Store only metadata (seed, K, phi, vol, cut) - NOT the full node list!
            TriageQueue.Push(ClusterInfo(phi, SeedNId, clusterSize, vol, cut));
            addedToQueue++;
          }
        }
      }
    }
  }
  
  printf("  Total scans: %d\n", totalScans);
  printf("  Clusters added to queue: %d\n", addedToQueue);
  printf("  Queue size: %d\n", TriageQueue.Size());
}

// Select target node based on Fiedler sign and PageRank
int SelectTarget(const PUNGraph& Graph, int SourceNId, 
                const TIntFltH& FiedlerVec, const TIntFltH& PageRankVec,
                const TIntSet& ClusterNodes) {
  // Validate source node
  if (!Graph->IsNode(SourceNId)) {
    return -1;  // Invalid source
  }

  double sourceFiedler = FiedlerVec.GetDat(SourceNId);
  int sourceSign = (sourceFiedler >= 0) ? 1 : -1;
  
  // Build candidate set (opposite Fiedler sign)
  TIntFltH Candidates;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    int nid = NI.GetId();
    if (ClusterNodes.IsKey(nid)) continue;  // Skip nodes in cluster
    if (nid == SourceNId) continue;  // Skip source itself

    double fiedler = FiedlerVec.GetDat(nid);
    int sign = (fiedler >= 0) ? 1 : -1;
    
    if (sign != sourceSign) {
      Candidates.AddDat(nid, PageRankVec.GetDat(nid));
    }
  }
  
  if (Candidates.Empty()) {
    // Fallback: any node not in cluster with high PageRank
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      int nid = NI.GetId();
      if (!ClusterNodes.IsKey(nid) && nid != SourceNId) {
        Candidates.AddDat(nid, PageRankVec.GetDat(nid));
      }
    }
  }
  
  if (Candidates.Empty()) {
    // Last resort: random node (excluding source)
    int target = Graph->GetRndNId();
    int attempts = 0;
    while (target == SourceNId && attempts < 100) {
      target = Graph->GetRndNId();
      attempts++;
    }
    return (target != SourceNId) ? target : -1;
  }
  
  // Weighted random selection
  double totalWeight = 0.0;
  for (int i = 0; i < Candidates.Len(); i++) {
    totalWeight += Candidates[i];
  }
  
  if (totalWeight <= 0.0) {
    // All weights are zero, use uniform random
    return Candidates.GetKey(TInt::Rnd.GetUniDevInt(Candidates.Len()));
  }

  double r = TFlt::Rnd.GetUniDev() * totalWeight;
  double cumWeight = 0.0;
  for (int i = 0; i < Candidates.Len(); i++) {
    cumWeight += Candidates[i];
    if (cumWeight >= r) {
      return Candidates.GetKey(i);
    }
  }
  
  return Candidates.GetKey(Candidates.Len() - 1);
}

// Phase III: Triage and Repair
void TriageAndRepair(PUNGraph& Graph,
                    ClusterPriorityQueue& TriageQueue,
                    const TIntFltH& FiedlerVec, const TIntFltH& PageRankVec,
                    int Budget, double PhiThreshold, TIntPrV& AddedEdges) {
  printf("\nPhase III: Triage and Repair\n");
  printf("  Budget: %d edges\n", Budget);
  printf("  Threshold: %.3f\n", PhiThreshold);
  
  int edgesAdded = 0;
  int recycled = 0;
  
  TLocClust Clust(Graph, 0.001);  // For recomputing clusters
  int skippedInvalidSeed = 0;
  int skippedInvalidTarget = 0;
  int skippedEmptyCluster = 0;

  while (Budget > 0 && !TriageQueue.IsEmpty()) {
    ClusterInfo cluster = TriageQueue.Pop();
    
    if (edgesAdded % 10 == 0 && edgesAdded > 0) {
      printf("  Progress: %d/%d edges added, Queue size: %d\n", 
             edgesAdded, edgesAdded + Budget, TriageQueue.Size());
    }
    
    // Validate seed node exists in graph
    if (!Graph->IsNode(cluster.SeedNId)) {
      skippedInvalidSeed++;
      continue;
    }

    // Recompute the cluster from seed node and K
    Clust.FindBestCut(cluster.SeedNId, cluster.K, 0.001);

    // Get the cluster nodes
    TIntV ClusterNIds;
    if (Clust.GetRndWalkSup() > 0) {
      // Find the cut with conductance closest to what we stored
      int bestIdx = 0;
      double minDiff = fabs(Clust.GetPhi(0) - cluster.Conductance);
      for (int i = 1; i < Clust.GetRndWalkSup(); i++) {
        double diff = fabs(Clust.GetPhi(i) - cluster.Conductance);
        if (diff < minDiff) {
          minDiff = diff;
          bestIdx = i;
        }
      }
      Clust.GetNIdV().GetSubValV(0, bestIdx, ClusterNIds);
    }

    if (ClusterNIds.Len() == 0) {
      skippedEmptyCluster++;
      continue;  // Skip if cluster couldn't be recomputed
    }

    // Source: seed node
    int source = cluster.SeedNId;
    
    // Build cluster node set for target selection
    TIntSet ClusterNodes(ClusterNIds.Len());
    for (int i = 0; i < ClusterNIds.Len(); i++) {
      ClusterNodes.AddKey(ClusterNIds[i]);
    }
    
    // Select target
    int target = SelectTarget(Graph, source, FiedlerVec, PageRankVec, ClusterNodes);
    
    // Validate target node
    if (target < 0 || !Graph->IsNode(target)) {
      skippedInvalidTarget++;
      continue;
    }

    // Skip self-loops
    if (source == target) {
      continue;
    }

    // Add edge if it doesn't exist
    if (!Graph->IsEdge(source, target)) {
      Graph->AddEdge(source, target);
      AddedEdges.Add(TIntPr(source, target));
      edgesAdded++;
      Budget--;
      
      // Recalculate conductance (approximate)
      int newVol = cluster.Volume + 1;
      int newCut = cluster.Cut + 1;  // Assuming target is outside cluster
      double newPhi = double(newCut) / double(newVol);
      
      // Recycle if still below threshold
      if (newPhi < PhiThreshold) {
        TriageQueue.Push(ClusterInfo(newPhi, source, cluster.K, newVol, newCut));
        recycled++;
      }
    }
  }
  
  printf("  Total edges added: %d\n", edgesAdded);
  printf("  Clusters recycled: %d\n", recycled);
  printf("  Final queue size: %d\n", TriageQueue.Size());
  if (skippedInvalidSeed > 0 || skippedInvalidTarget > 0 || skippedEmptyCluster > 0) {
    printf("  Skipped clusters:\n");
    if (skippedInvalidSeed > 0) {
      printf("    Invalid seed nodes: %d\n", skippedInvalidSeed);
    }
    if (skippedEmptyCluster > 0) {
      printf("    Empty clusters: %d\n", skippedEmptyCluster);
    }
    if (skippedInvalidTarget > 0) {
      printf("    Invalid targets: %d\n", skippedInvalidTarget);
    }
  }
}

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("SNAP-Guided Spectral Repair. build: %s, %s. Time: %s", 
                         __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  
  // Command line arguments
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", 
                                           "Input undirected graph (one edge per line)");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "repaired", 
                                            "Output file prefix");
  const int Budget = Env.GetIfArgPrefixInt("-b:", 100, 
                                           "Edge budget (number of edges to add)");
  const double PhiThreshold = Env.GetIfArgPrefixFlt("-t:", 0.1,
                                                     "Conductance threshold for valley detection");
  const int KMin = Env.GetIfArgPrefixInt("-kmin:", 10, 
                                        "Minimum cluster size K");
  const int KMax = Env.GetIfArgPrefixInt("-kmax:", 100000,
                                        "Maximum cluster size K");
  const double KFac = Env.GetIfArgPrefixFlt("-kfac:", 1.2, 
                                           "K growth factor");
  const int Coverage = Env.GetIfArgPrefixInt("-c:", 10, 
                                            "Coverage (scans per node)");
  const bool RunNCP = Env.GetIfArgPrefixBool("-ncp:", true, 
                                             "Run NCP plot after repair");
  
  // New options for queue management
  const bool BuildQueue = Env.GetIfArgPrefixBool("-buildq:", true,
                                                 "Build priority queue (Phase I & II)");
  const bool AddEdges = Env.GetIfArgPrefixBool("-repair:", true,
                                               "Add edges (Phase III)");
  const TStr SaveQueueFNm = Env.GetIfArgPrefixStr("-saveq:", "",
                                                  "Save priority queue to file (empty = don't save)");
  const TStr LoadQueueFNm = Env.GetIfArgPrefixStr("-loadq:", "",
                                                  "Load priority queue from file (empty = don't load)");

  printf("SNAP-Guided Spectral Repair\n");
  printf("============================\n");
  printf("Input graph: %s\n", InFNm.CStr());
  printf("Output prefix: %s\n", OutFNm.CStr());
  printf("Budget: %d edges\n", Budget);
  printf("Phi threshold: %.3f\n", PhiThreshold);
  printf("Build Queue: %s\n", BuildQueue ? "Yes" : "No");
  printf("Add Edges: %s\n", AddEdges ? "Yes" : "No");
  if (!SaveQueueFNm.Empty()) {
    printf("Save Queue to: %s\n", SaveQueueFNm.CStr());
  }
  if (!LoadQueueFNm.Empty()) {
    printf("Load Queue from: %s\n", LoadQueueFNm.CStr());
  }
  printf("\n");
  
  // Load graph
  printf("Loading graph...\n");
  PUNGraph Graph = TSnap::GetMxWcc(TSnap::LoadEdgeList<PUNGraph>(InFNm, 0, 1));
  printf("  Nodes: %d, Edges: %d\n", Graph->GetNodes(), Graph->GetEdges());
  
  // Initialize priority queue
  ClusterPriorityQueue TriageQueue;

  // Phase I & II: Build Queue (or load from file)
  if (!LoadQueueFNm.Empty()) {
    printf("\nLoading priority queue from file...\n");
    TriageQueue.Load(LoadQueueFNm);
    printf("  Loaded queue with %d clusters\n", TriageQueue.Size());
  } else if (BuildQueue) {
    // Phase I: Global Mapping
    printf("\nPhase I: Global Mapping\n");
    TIntFltH FiedlerVec, PageRankVec;
    ComputeFiedlerVector(Graph, FiedlerVec);
    ComputeGlobalPageRank(Graph, PageRankVec);

    // Phase II: Build Triage Queue
    BuildTriageQueue(Graph, TriageQueue, KMin, KMax, KFac, Coverage, PhiThreshold);
  }

  // Save queue if requested
  if (!SaveQueueFNm.Empty() && TriageQueue.Size() > 0) {
    printf("\nSaving priority queue to file...\n");
    TriageQueue.Save(SaveQueueFNm);
    printf("  Queue saved to: %s\n", SaveQueueFNm.CStr());
    printf("  Queue size: %d clusters\n", TriageQueue.Size());
  }

  // Phase III: Triage and Repair (if requested)
  TIntPrV AddedEdges;
  if (AddEdges && TriageQueue.Size() > 0) {
    // Need to recompute Fiedler and PageRank for target selection
    printf("\nComputing global vectors for target selection...\n");
    TIntFltH FiedlerVec, PageRankVec;
    ComputeFiedlerVector(Graph, FiedlerVec);
    ComputeGlobalPageRank(Graph, PageRankVec);

    TriageAndRepair(Graph, TriageQueue, FiedlerVec, PageRankVec, Budget, PhiThreshold, AddedEdges);
  }

  // Save results (if edges were added)
  if (AddedEdges.Len() > 0) {
    printf("\nSaving results...\n");

    // Save added edges
    TStr EdgesFNm = OutFNm + "_edges.txt";
    FILE* FEdges = fopen(EdgesFNm.CStr(), "wt");
    fprintf(FEdges, "# Added edges from SNAP-Guided Spectral Repair\n");
    fprintf(FEdges, "# Format: source target\n");
    fprintf(FEdges, "# Total edges added: %d\n", AddedEdges.Len());
    for (int i = 0; i < AddedEdges.Len(); i++) {
      fprintf(FEdges, "%d\t%d\n", AddedEdges[i].Val1.Val, AddedEdges[i].Val2.Val);
    }
    fclose(FEdges);
    printf("  Edges saved to: %s\n", EdgesFNm.CStr());

    // Save modified graph
    TStr GraphFNm = OutFNm + "_graph.txt";
    TSnap::SaveEdgeList(Graph, GraphFNm.CStr());
    printf("  Modified graph saved to: %s\n", GraphFNm.CStr());
    printf("  Final graph: Nodes=%d, Edges=%d\n", Graph->GetNodes(), Graph->GetEdges());

    // Optionally run NCP plot
    if (RunNCP) {
      printf("\nRunning NCP plot on modified graph...\n");
      TLocClust::PlotNCP(Graph, OutFNm + "_ncp",
                        "SNAP-Guided Spectral Repair Results",
                        false, false, 10, Mega(100), 10, false, false, -1);
    }
  }
  
  Catch
  printf("\nTotal run time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  
  return 0;
}
