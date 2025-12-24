#ifndef snap_cluster_manager_hpp
#define snap_cluster_manager_hpp

#include "stdafx.h"

// Enhanced ClusterInfo structure for dynamic repair
// Tracks state needed for local re-computation and self-cleaning
struct ClusterInfo {
    double Conductance;
    int SeedNId;
    int K;              // Original size parameter for local sweep
    int LastUpdated;    // Total edges added to graph when this was last checked

    ClusterInfo() : Conductance(1.0), SeedNId(-1), K(0), LastUpdated(0) {}
    
    ClusterInfo(double phi, int seed, int k, int time) 
        : Conductance(phi), SeedNId(seed), K(k), LastUpdated(time) {}

    // Priority Queue: Lowest conductance (worst bottlenecks) at the top
    bool operator<(const ClusterInfo& other) const {
        return Conductance > other.Conductance; 
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

#endif  // snap_cluster_manager_hpp
