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

// Compute Fiedler vector (second smallest eigenvector of Laplacian)
void ComputeFiedlerVector(const PUNGraph& Graph, TIntFltH& FiedlerVec) {
  printf("Computing Fiedler vector...\n");
  // For now, use a simple approximation based on degree centrality
  // A full implementation would use Lanczos or power iteration
  FiedlerVec.Clr();
  
  // Simple heuristic: assign positive values to high-degree nodes
  // and negative values to low-degree nodes (this is a placeholder)
  TIntV DegV;
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    DegV.Add(NI.GetDeg());
  }
  DegV.Sort();
  double median = (DegV.Len() > 0) ? double(DegV[DegV.Len()/2].Val) : 1.0;
  
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    int nid = NI.GetId();
    double val = (NI.GetDeg() - median) / (median + 1.0);
    FiedlerVec.AddDat(nid, val);
  }
  printf("  Fiedler vector computed for %d nodes\n", FiedlerVec.Len());
}

// Compute Global PageRank
void ComputeGlobalPageRank(const PUNGraph& Graph, TIntFltH& PageRankVec, double Alpha=0.85) {
  printf("Computing Global PageRank...\n");
  const int MaxIter = 100;
  const double Eps = 1e-6;
  const int N = Graph->GetNodes();
  
  PageRankVec.Clr();
  TIntFltH OldPR;
  
  // Initialize
  for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    PageRankVec.AddDat(NI.GetId(), 1.0 / N);
  }
  
  // Power iteration
  for (int iter = 0; iter < MaxIter; iter++) {
    OldPR = PageRankVec;
    double diff = 0.0;
    
    for (TUNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      int nid = NI.GetId();
      double sum = 0.0;
      
      // Sum contributions from neighbors
      for (int e = 0; e < NI.GetInDeg(); e++) {
        int nbr = NI.GetInNId(e);
        int nbrDeg = Graph->GetNI(nbr).GetOutDeg();
        if (nbrDeg > 0) {
          sum += OldPR.GetDat(nbr) / nbrDeg;
        }
      }
      
      double newPR = (1.0 - Alpha) / N + Alpha * sum;
      PageRankVec.AddDat(nid, newPR);
      diff += fabs(newPR - OldPR.GetDat(nid));
    }
    
    if (diff < Eps) {
      printf("  PageRank converged after %d iterations\n", iter+1);
      break;
    }
  }
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
                     TVec<ClusterInfo>& TriageVec,
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
            
            TriageVec.Add(ClusterInfo(phi, SeedNId, vol, cut, ClusterNIds));
            addedToQueue++;
          }
        }
      }
    }
  }
  
  printf("  Total scans: %d\n", totalScans);
  printf("  Clusters added to queue: %d\n", addedToQueue);
  
  // Sort by conductance (ascending)
  TriageVec.Sort(true);
  printf("  Queue size: %d\n", TriageVec.Len());
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
                    TVec<ClusterInfo>& TriageVec,
                    const TIntFltH& FiedlerVec, const TIntFltH& PageRankVec,
                    int Budget, double PhiThreshold, TIntPrV& AddedEdges) {
  printf("\nPhase III: Triage and Repair\n");
  printf("  Budget: %d edges\n", Budget);
  printf("  Threshold: %.3f\n", PhiThreshold);
  
  int edgesAdded = 0;
  int recycled = 0;
  int idx = 0;
  
  while (Budget > 0 && idx < TriageVec.Len()) {
    ClusterInfo cluster = TriageVec[idx];
    idx++;
    
    if (edgesAdded % 10 == 0 && edgesAdded > 0) {
      printf("  Progress: %d/%d edges added, Queue remaining: %d\n", 
             edgesAdded, edgesAdded + Budget, TriageVec.Len() - idx);
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
      
      // Recycle if still below threshold - add back to queue
      if (newPhi < PhiThreshold) {
        TriageVec.Add(ClusterInfo(newPhi, source, newVol, newCut, cluster.NodeIds));
        recycled++;
      }
    }
  }
  
  printf("  Total edges added: %d\n", edgesAdded);
  printf("  Clusters recycled: %d\n", recycled);
  printf("  Final queue size: %d\n", TriageVec.Len() - idx);
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
  TVec<ClusterInfo> TriageVec;
  BuildTriageQueue(Graph, TriageVec, KMin, KMax, KFac, Coverage, PhiThreshold);
  
  // Phase III: Triage and Repair
  TIntPrV AddedEdges;
  TriageAndRepair(Graph, TriageVec, FiedlerVec, PageRankVec, Budget, PhiThreshold, AddedEdges);
  
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
