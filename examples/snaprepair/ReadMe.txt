# SNAP-Guided Spectral Repair

This tool implements the SNAP-Guided Spectral Repair algorithm for systematically flattening the Network Community Profile (NCP) curve through targeted graph rewiring.

## Overview

The SNAP-Guided Spectral Repair algorithm addresses the problem of "local traps" or low-conductance clusters in networks by:
1. Identifying clusters with poor connectivity (low conductance)
2. Adding edges strategically to improve global mixing properties
3. Using spectral and PageRank guidance to select optimal rewiring targets

The result is a network with a flatter NCP curve, indicating more uniform community structure across scales.

## Algorithm Overview

The algorithm operates in three phases:

### Phase I: Global Mapping
Computes global network properties to guide rewiring:
- **Fiedler Vector**: Second eigenvector of the graph Laplacian, provides spectral coordinates
- **Global PageRank**: Identifies important "hub" nodes in the network

### Phase II: Robust SNAP Scanner
Uses the NCP algorithm to identify problematic clusters:
- Performs multiple Approximate PageRank (APR) runs at various scales
- Generates conductance profiles via sweep cuts
- Detects "valleys" (clusters with conductance below threshold)
- Populates a priority queue ordered by conductance (worst first)
- Deduplicates clusters to avoid redundant processing

### Phase III: Triage and Repair
Systematically improves the network:
- Processes clusters in order of increasing conductance
- For each cluster, selects source node (seed) and target node
- Target selection uses:
  - **Spectral Opposition**: Chooses nodes on opposite side of Fiedler cut
  - **PageRank Weighting**: Prefers high-PageRank hubs for better mixing
  - **Randomization**: Maintains diversity to avoid star topology
- Adds edge and recalculates cluster conductance
- Recycles cluster into queue if still below threshold

## Building

```bash
cd examples/snaprepair
make
```

## Usage

```bash
./snaprepair -i:<input_graph> -o:<output_prefix> [options]
```

### Required Arguments

- `-i:<file>` - Input undirected graph file (one edge per line, tab/space separated)

### Optional Arguments

- `-o:<prefix>` - Output file prefix (default: "repaired")
- `-b:<int>` - Edge budget, number of edges to add (default: 100)
- `-t:<float>` - Conductance threshold for valley detection (default: 0.2)
- `-kmin:<int>` - Minimum cluster size K (default: 10)
- `-kmax:<int>` - Maximum cluster size K (default: 100)
- `-kfac:<float>` - K growth factor (default: 1.2)
- `-c:<int>` - Coverage, scans per node (default: 10)
- `-ncp:<bool>` - Run NCP plot after repair (default: true)

## Output Files

The algorithm produces:

1. `<prefix>_edges.txt` - List of added edges (one per line: source target)
2. `<prefix>_graph.txt` - Complete modified graph edge list
3. `<prefix>_ncp.*` - NCP plot files (if -ncp:T)

## Example

Basic usage:
```bash
# Run with 200 edge budget on a graph
./snaprepair -i:../as20graph.txt -o:repaired_as20 -b:200 -t:0.2

# Then compare with original using ncpplot
cd ../ncpplot
./ncpplot -i:../../as20graph.txt -o:original
./ncpplot -i:../snaprepair/repaired_as20_graph.txt -o:repaired
```

Complete demo workflow:
```bash
# Run the demo script
./demo.sh

# This will:
# 1. Run the repair algorithm
# 2. Generate NCP plots for both original and repaired graphs
# 3. Show comparison statistics
```

Quick test with small parameters:
```bash
./snaprepair -i:../as20graph.txt -o:test -b:50 -kmin:10 -kmax:50 -ncp:F
```

## Algorithm Parameters

- **Budget (-b)**: Number of edges to add. Higher values = more aggressive flattening
- **Phi Threshold (-t)**: Conductance threshold for identifying "valleys". Lower = stricter
- **KMin, KMax (-kmin, -kmax)**: Range of cluster sizes to scan. Should cover the range where dips appear in the NCP plot
- **KFac (-kfac)**: Growth factor for K. Values around 1.2-1.5 provide good coverage
- **Coverage (-c)**: How many times to scan each potential cluster size. Higher = more thorough but slower

## Notes

- The algorithm works on the largest weakly connected component of the input graph
- Added edges connect nodes on opposite sides of the Fiedler cut (spectrally distant)
- Target nodes are selected with probability proportional to their PageRank (preferring hubs)
- Clusters can be "recycled" back into the priority queue if they remain problematic after rewiring
- The Fiedler vector computation uses a degree-based heuristic that provides reasonable spectral partitioning without requiring external eigenvalue libraries

## Expected Results

After running the algorithm, you should observe:

1. **Flatter NCP Curve**: The conductance at small cluster sizes (k=10-100) should increase
2. **Reduced Valleys**: Sharp dips in the original NCP plot should be smoothed out
3. **Better Mixing**: Random walks spread more evenly throughout the network
4. **Maintained Structure**: Large-scale structure is preserved while local connectivity improves

## Algorithm Parameters and Tuning

**Budget (-b)**: 
- Typical values: 50-500 for small graphs (< 10k nodes), 100-1000 for medium graphs
- Higher values = more aggressive flattening
- Rule of thumb: Start with 1-5% of original edge count

**Phi Threshold (-t)**: 
- Default: 0.2 is generally effective
- Lower values (0.1-0.15): More selective, targets only worst clusters
- Higher values (0.25-0.3): More aggressive, targets broader range of clusters

**K Range (-kmin, -kmax)**: 
- Should cover the range where dips appear in the NCP plot
- Default: 10-100 works well for most networks
- For large networks, extend -kmax to 1000 or more
- Use -kfac to control sampling density (1.2-1.5 recommended)

**Coverage (-c)**: 
- Controls thoroughness of scanning
- Default: 10 provides good balance
- Higher values (15-20): More comprehensive but slower
- Lower values (5-8): Faster but may miss some clusters

## Troubleshooting

**No clusters found**: 
- Your graph may already have good mixing properties
- Try increasing -t (phi threshold) to be more inclusive
- Check that -kmin and -kmax cover the right range

**Algorithm runs but no improvement visible**:
- Increase budget (-b)
- Lower threshold (-t) to target worse clusters
- Ensure NCP plots are generated with same parameters

**Very slow execution**:
- Reduce -kmax to limit cluster sizes scanned
- Reduce -c (coverage) for faster but less thorough scanning
- Use smaller budget for testing

## Technical Details

**Complexity**:
- Phase I (PageRank): O(|E| * iterations) ≈ O(|E|)
- Phase II (Scanning): O(Coverage * |E| * log(Kmax/Kmin))
- Phase III (Repair): O(Budget * (|V| + |E|))
- Overall: O(Coverage * |E| * log(K) + Budget * |V|)

**Memory Usage**:
- Priority queue: O(number of detected clusters)
- Graph storage: O(|V| + |E|)
- Typical: < 1GB for graphs with 100k nodes

**Randomization**:
- The algorithm uses randomization in:
  1. Seed node selection during scanning
  2. Target node selection (weighted by PageRank)
- Results may vary slightly between runs
- For reproducibility, set TRnd seed before running
