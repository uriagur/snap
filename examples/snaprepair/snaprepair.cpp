#include "stdafx.h"
#include "ncp.h"

// Structure to hold cluster information for the priority queue
struct ClusterInfo {
  double Conductance;
  int SeedNId;
  int Volume;
  int Cut;
  TIntV NodeIds;  // Node IDs in the cluster
  
  ClusterInfo() : Conductance(1.0), SeedNId(-1), Volume(0), Cut(0) {}
  
  ClusterInfo(double phi, int seed, int vol, int cut, const TIntV& nids) :
    Conductance(phi), SeedNId(seed), Volume(vol), Cut(cut), NodeIds(nids) {}
  
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
};


// Compute Fiedler vector approximation using degree-based heuristic
// The Fiedler vector (second eigenvector of graph Laplacian) provides spectral coordinates.
// This approximation uses degree centrality normalized around the median degree.
// Nodes with above-median degree get positive values, below-median get negative values.
// This provides a reasonable spectral partition for target selection.
void ComputeFiedlerVector(const PUNGraph& Graph, TIntFltH& FiedlerVec) {
  printf("Computing Fiedler vector approximation...\n");
  FiedlerVec.Clr();
  
  // Collect all node degrees
  TIntV DegV;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    DegV.Add(NI.GetDeg());
  }
  
  // Sort to find median
  DegV.Sort();
  double median = (DegV.Len() > 0) ? double(DegV[DegV.Len()/2].Val) : 1.0;
  
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
          // Extract cluster
          TIntV ClusterNIds;
          Clust.GetNIdV().GetSubValV(0, i, ClusterNIds);
          
          // Check if we've seen this cluster
          TUInt64 hash = HashNodeSet(ClusterNIds);
          if (!SeenClusters.IsKey(hash)) {
            SeenClusters.AddDat(hash, true);
            
            int vol = Clust.GetVol(i);
            int cut = Clust.GetCut(i);
            
            TriageQueue.Push(ClusterInfo(phi, SeedNId, vol, cut, ClusterNIds));
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
  double sourceFiedler = FiedlerVec.GetDat(SourceNId);
  int sourceSign = (sourceFiedler >= 0) ? 1 : -1;
  
  // Build candidate set (opposite Fiedler sign)
  TIntFltH Candidates;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    int nid = NI.GetId();
    if (ClusterNodes.IsKey(nid)) continue;  // Skip nodes in cluster
    
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
      if (!ClusterNodes.IsKey(nid)) {
        Candidates.AddDat(nid, PageRankVec.GetDat(nid));
      }
    }
  }
  
  if (Candidates.Empty()) {
    return Graph->GetRndNId();
  }
  
  // Weighted random selection
  double totalWeight = 0.0;
  for (int i = 0; i < Candidates.Len(); i++) {
    totalWeight += Candidates[i];
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
  
  while (Budget > 0 && !TriageQueue.IsEmpty()) {
    ClusterInfo cluster = TriageQueue.Pop();
    
    if (edgesAdded % 10 == 0 && edgesAdded > 0) {
      printf("  Progress: %d/%d edges added, Queue size: %d\n", 
             edgesAdded, edgesAdded + Budget, TriageQueue.Size());
    }
    
    // Source: seed node
    int source = cluster.SeedNId;
    
    // Build cluster node set for target selection
    TIntSet ClusterNodes(cluster.NodeIds.Len());
    for (int i = 0; i < cluster.NodeIds.Len(); i++) {
      ClusterNodes.AddKey(cluster.NodeIds[i]);
    }
    
    // Select target
    int target = SelectTarget(Graph, source, FiedlerVec, PageRankVec, ClusterNodes);
    
    // Add edge if it doesn't exist
    if (!Graph->IsEdge(source, target)) {
      Graph->AddEdge(source, target);
      AddedEdges.Add(TIntPr(source, target));
      edgesAdded++;
      Budget--;
      
      // Recalculate conductance
      int newVol = cluster.Volume + 1;
      int newCut = cluster.Cut + 1;  // Assuming target is outside cluster
      double newPhi = double(newCut) / double(newVol);
      
      // Recycle if still below threshold - push back to priority queue
      if (newPhi < PhiThreshold) {
        TriageQueue.Push(ClusterInfo(newPhi, source, newVol, newCut, cluster.NodeIds));
        recycled++;
      }
    }
  }
  
  printf("  Total edges added: %d\n", edgesAdded);
  printf("  Clusters recycled: %d\n", recycled);
  printf("  Final queue size: %d\n", TriageQueue.Size());
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
  const bool RunNCP = Env.GetIfArgPrefixBool("-ncp:", true, 
                                             "Run NCP plot after repair");
  
  printf("SNAP-Guided Spectral Repair\n");
  printf("============================\n");
  printf("Input graph: %s\n", InFNm.CStr());
  printf("Output prefix: %s\n", OutFNm.CStr());
  printf("Budget: %d edges\n", Budget);
  printf("Phi threshold: %.3f\n", PhiThreshold);
  printf("\n");
  
  // Load graph
  printf("Loading graph...\n");
  PUNGraph Graph = TSnap::GetMxWcc(TSnap::LoadEdgeList<PUNGraph>(InFNm, 0, 1));
  printf("  Nodes: %d, Edges: %d\n", Graph->GetNodes(), Graph->GetEdges());
  
  // Phase I: Global Mapping
  printf("\nPhase I: Global Mapping\n");
  TIntFltH FiedlerVec, PageRankVec;
  ComputeFiedlerVector(Graph, FiedlerVec);
  ComputeGlobalPageRank(Graph, PageRankVec);
  
  // Phase II: Build Triage Queue
  ClusterPriorityQueue TriageQueue;
  BuildTriageQueue(Graph, TriageQueue, KMin, KMax, KFac, Coverage, PhiThreshold);
  
  // Phase III: Triage and Repair
  TIntPrV AddedEdges;
  TriageAndRepair(Graph, TriageQueue, FiedlerVec, PageRankVec, Budget, PhiThreshold, AddedEdges);
  
  // Save results
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
  
  Catch
  printf("\nTotal run time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  
  return 0;
}
