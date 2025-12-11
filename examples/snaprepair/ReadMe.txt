# SNAP-Guided Spectral Repair

This tool implements the SNAP-Guided Spectral Repair algorithm for systematically flattening the Network Community Profile (NCP) curve through targeted graph rewiring.

## Algorithm Overview

The algorithm operates in three phases:

1. **Phase I: Global Mapping** - Computes Fiedler vector and Global PageRank to establish a coordinate system for rewiring
2. **Phase II: Robust SNAP Scanner** - Uses NCP scanning to detect low-conductance clusters (local traps) and populates a priority queue
3. **Phase III: Triage and Repair** - Adds edges between spectrally distant, high-PageRank nodes to flatten the NCP curve

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

```bash
# Run with 200 edge budget on a graph
./snaprepair -i:../as20graph.txt -o:repaired_as20 -b:200 -t:0.2

# Then compare with original using ncpplot
cd ../ncpplot
./ncpplot -i:../../as20graph.txt -o:original
./ncpplot -i:../snaprepair/repaired_as20_graph.txt -o:repaired
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
