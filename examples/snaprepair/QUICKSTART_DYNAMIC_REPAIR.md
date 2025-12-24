# Quick Start Guide - Dynamic Repair Algorithm

## Building the Algorithm

### Step 1: Navigate to the snaprepair directory
```bash
cd examples/snaprepair
```

### Step 2: Build the dynamic repair test program
```bash
# Compile the test program
g++ -std=c++98 -Wall -O3 -DNDEBUG -fopenmp \
    -I../../glib-core -I../../snap-core -I../../snap-adv \
    test_dynamic_repair.cpp RepairStrategies.cpp DynamicRepair.cpp \
    ../../snap-adv/ncp.cpp -o test_dynamic_repair \
    ../../snap-core/Snap.o -lrt
```

Or simply use the existing makefile (which builds snaprepair):
```bash
make
```

## Running the Algorithm

### Basic Usage
```bash
./test_dynamic_repair -i:../as20graph.txt -b:100 -t:0.2
```

### Command-line Parameters
- `-i:` Input graph file (edge list format, one edge per line)
- `-o:` Output file prefix (default: dynamic_repair)
- `-b:` Edge budget - number of edges to add (default: 100)
- `-t:` Conductance threshold for valley detection (default: 0.2)
- `-kmin:` Minimum cluster size K (default: 10)
- `-kmax:` Maximum cluster size K (default: 100)
- `-kfac:` K growth factor (default: 1.2)
- `-c:` Coverage - scans per node (default: 10)

### Example Runs

**Quick Test (small budget):**
```bash
./test_dynamic_repair -i:../as20graph.txt -b:10 -t:0.2 -kmin:10 -kmax:50 -c:2
```

**Standard Run:**
```bash
./test_dynamic_repair -i:../as20graph.txt -b:100 -t:0.2
```

**Aggressive Repair (more edges, lower threshold):**
```bash
./test_dynamic_repair -i:../as20graph.txt -b:200 -t:0.15 -kmax:200
```

## How It Works

The algorithm follows a two-phase approach as designed:

### Phase 1: Build Priority Queue
The algorithm first scans the graph to identify problematic clusters:

```
Building Triage Queue for Dynamic Repair
  KMin=10, KMax=100, KFac=1.20, Coverage=10, PhiThreshold=0.200
  K=10, Eps=0.010000, Runs=2780
  K=13, Eps=0.007692, Runs=2138
  ...
  Total scans: 64524
  Clusters added to queue: 134
  Queue size: 134
```

**What happens:**
1. Loops through cluster sizes K from KMin to KMax (growing by KFac)
2. For each K, runs multiple local clustering operations (Runs determined by Coverage)
3. Uses `TLocClust::FindBestCut()` to compute conductance
4. Adds clusters with conductance < threshold to priority queue
5. Queue orders by conductance (worst bottlenecks first)

### Phase 2: Dynamic Repair Loop
Then iterates through the priority queue to fix clusters:

```
Running Dynamic Repair
  Budget: 100 edges
  Threshold: 0.200
  Progress: 10 edges added, 20 clusters processed, Queue size: 114
  ...
  Total edges added: 100
  Clusters processed: 133
  Clusters fixed: 123
  Clusters recycled: 10
```

**What happens:**
1. **Pop** worst cluster from priority queue
2. **VERIFY**: Re-run local sweep to check if still bad
3. **REPAIR**: If still bad (φ < threshold), add edge using selected strategy
4. **RE-EVALUATE**: Immediately check new conductance after edge addition
5. **RECYCLE**: If still bad, push back to queue with updated conductance
6. **SELF-CLEAN**: If fixed (φ ≥ threshold), discard cluster
7. Repeat until budget exhausted or queue empty

## Understanding the Output

### Console Output
```
Dynamic Repair Test
===================
Input graph: ../as20graph.txt
Output prefix: dynamic_repair
Budget: 100 edges
Phi threshold: 0.200

Loading graph...
  Nodes: 6474, Edges: 13895

Building Triage Queue for Dynamic Repair
  [Phase 1 output...]
  
Running Dynamic Repair
  [Phase 2 output...]
  
Saving results...
  Modified graph saved to: dynamic_repair_graph.txt
  Final graph: Nodes=6474, Edges=13995

Total run time: 0.03s
```

### Output Files

**`dynamic_repair_graph.txt`**
- The modified graph with added edges
- Format: same as input (edge list)
- Can be used as input to other SNAP tools

### Key Metrics

**Clusters Fixed**: Clusters that had conductance improved to ≥ threshold
**Clusters Recycled**: Clusters that needed multiple edges (pushed back to queue)
**Self-Cleaning Rate**: (Fixed / Processed) - typically 95%+

## Architecture Overview

The implementation uses three modular components:

### 1. ClusterManager.hpp
- `ClusterInfo` struct: Stores conductance, seed node, K, and timestamp
- `ClusterPriorityQueue`: Min-heap ordering by conductance

### 2. RepairStrategies.cpp/h
- `SelectTarget()`: Pluggable target selection
- `STRATEGY_RANDOM`: Boundary-biased random (currently used)
- Extensible for other strategies (e.g., spectral/Fiedler)

### 3. DynamicRepair.cpp/h
- `RunDynamicRepair()`: Main loop with 3-step process
- Local verification, edge addition, immediate re-evaluation
- Self-cleaning: discards fixed clusters automatically

## Integration with Your Code

To use the dynamic repair algorithm in your own code:

```cpp
#include "ClusterManager.hpp"
#include "RepairStrategies.h"
#include "DynamicRepair.h"

// Load your graph
PUNGraph Graph = TSnap::LoadEdgeList<PUNGraph>("mygraph.txt", 0, 1);

// Phase 1: Build priority queue
ClusterPriorityQueue TriageQueue;
BuildTriageQueueForDynamicRepair(Graph, TriageQueue, 
                                 10, 100, 1.2, 10, 0.2);

// Phase 2: Run dynamic repair
RunDynamicRepair(Graph, TriageQueue, 100, 0.2);

// Save results
TSnap::SaveEdgeList(Graph, "repaired_graph.txt");
```

## Performance Tips

**For faster execution:**
- Reduce Coverage (`-c:`) - fewer scans per K
- Reduce KMax (`-kmax:`) - scan smaller clusters only
- Use smaller budget for testing

**For more thorough repair:**
- Increase Coverage (`-c:`) - more comprehensive scanning
- Increase KMax (`-kmax:`) - include larger clusters
- Lower threshold (`-t:`) - target more clusters

## Troubleshooting

**No clusters added to queue:**
- Graph may already have good conductance
- Try increasing threshold (`-t:`)
- Check KMin/KMax range is appropriate for your graph

**Too many clusters recycled:**
- Normal if clusters are very problematic
- Consider increasing budget (`-b:`)
- May need multiple passes

**Algorithm is slow:**
- Reduce Coverage (`-c:`)
- Reduce KMax (`-kmax:`)
- Use smaller budget for testing

## Comparison with Original snaprepair

| Feature | Original snaprepair | Dynamic Repair |
|---------|-------------------|----------------|
| Phase 1 | Global mapping (Fiedler, PageRank) | Local clustering scan |
| Phase 2 | Build triage queue | Build priority queue |
| Phase 3 | Static repair | Dynamic repair with recycling |
| Verification | None (one-shot) | Per-cluster re-verification |
| Recycling | Simple (estimate) | Immediate re-evaluation |
| Self-cleaning | No | Yes (auto-discard fixed clusters) |
| Target selection | Spectral/PageRank | Pluggable strategies |

Both algorithms work well - choose based on your needs:
- **snaprepair**: Global spectral information, good for overall structure
- **dynamic repair**: Local refinement, efficient recycling, modular

## Next Steps

1. **Run the test**: Try the examples above
2. **Read the docs**: See `DYNAMIC_REPAIR.md` for architecture details
3. **Customize**: Implement your own repair strategy in RepairStrategies.cpp
4. **Integrate**: Use the components in your own code

For more details, see:
- `DYNAMIC_REPAIR.md` - Full documentation
- `IMPLEMENTATION_SUMMARY.md` - Technical reference
- `test_dynamic_repair.cpp` - Complete example
