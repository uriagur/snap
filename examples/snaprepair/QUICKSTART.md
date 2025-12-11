# Quick Start Guide - SNAP-Guided Spectral Repair

## 1-Minute Quick Start

```bash
# Navigate to the snaprepair directory
cd examples/snaprepair

# Build the tool (if not already built)
make

# Run the demo to see everything in action
./demo.sh
```

That's it! The demo will:
- Run the repair algorithm on the AS20 graph
- Add 100 edges strategically
- Generate NCP plots for comparison
- Show you the results

## 5-Minute Tutorial

### Step 1: Run the algorithm on your graph
```bash
./snaprepair -i:your_graph.txt -o:my_repair -b:100
```

This creates:
- `my_repair_edges.txt` - The 100 edges that were added
- `my_repair_graph.txt` - Your modified graph

### Step 2: Compare NCP plots
```bash
./compare_ncp.sh your_graph.txt my_repair
```

This generates NCP plots for both original and repaired graphs.

### Step 3: Examine the results
- Check `my_repair_edges.txt` to see which edges were added
- The algorithm targets low-conductance clusters
- Added edges connect spectrally distant, important nodes

## What Does It Do?

The algorithm finds "problematic" parts of your network (low-conductance clusters that trap random walks) and adds edges to improve mixing:

1. **Scans** your graph at multiple scales to find clusters
2. **Identifies** clusters with poor conductance (φ < 0.2 by default)
3. **Repairs** by adding edges to spectrally distant, high-PageRank nodes
4. **Result**: Flatter NCP curve = better global mixing

## Common Use Cases

### Use Case 1: Quick Test
```bash
./snaprepair -i:../as20graph.txt -o:test -b:50 -kmin:10 -kmax:50 -ncp:F
```
Fast test with small parameters.

### Use Case 2: Aggressive Repair
```bash
./snaprepair -i:my_graph.txt -o:aggressive -b:500 -t:0.25
```
Adds many edges, targets broader range of clusters.

### Use Case 3: Conservative Repair
```bash
./snaprepair -i:my_graph.txt -o:conservative -b:100 -t:0.15
```
Fewer edges, only targets worst clusters.

## Parameter Cheat Sheet

| Want to... | Adjust this parameter |
|------------|----------------------|
| Add more edges | Increase `-b:` (budget) |
| Target more clusters | Increase `-t:` (threshold) |
| Target only worst clusters | Decrease `-t:` (threshold) |
| Scan larger clusters | Increase `-kmax:` |
| Make it faster | Reduce `-c:` (coverage) |
| Make it more thorough | Increase `-c:` (coverage) |

## Troubleshooting

**Q: The algorithm runs but I don't see any improvement**
- Try increasing the budget (`-b:`)
- Lower the threshold (`-t:`) to target worse clusters
- Make sure your graph actually has low-conductance regions (check original NCP plot)

**Q: It's taking too long**
- Reduce `-kmax:` to scan smaller cluster sizes only
- Reduce `-c:` for less thorough scanning
- Use smaller budget for testing

**Q: No clusters were found**
- Your graph may already have good mixing
- Try increasing `-t:` to be more inclusive
- Check that `-kmin:` and `-kmax:` cover the right range

## Understanding the Output

When you run the algorithm, you'll see:

```
Phase I: Global Mapping
  Computing Fiedler vector...     # Spectral coordinates
  Computing Global PageRank...    # Hub importance

Phase II: Building Triage Queue
  K=10, Eps=0.01, Runs=13892     # Scanning at scale K=10
  K=13, Eps=0.008, Runs=10682    # Scanning at scale K=13
  ...
  Clusters added to queue: 134    # Found 134 problems

Phase III: Triage and Repair
  Progress: 10/100 edges added   # Adding edges
  Total edges added: 100         # Done
  Clusters recycled: 77          # Some needed multiple edges
```

## Files Generated

After running, you get:

1. **`*_edges.txt`**
   - List of edges added
   - Format: `source_node\ttarget_node`
   - One edge per line

2. **`*_graph.txt`**
   - Complete modified graph
   - Same format as input
   - Can be used as input to other SNAP tools

3. **NCP plot files** (if `-ncp:T`)
   - `.tab` files with numerical data
   - `.plt` files for gnuplot
   - `.png` files (if gnuplot available)

## Getting Help

- **Full documentation**: See `ReadMe.txt`
- **Technical details**: See `IMPLEMENTATION.md`
- **Run the demo**: `./demo.sh`
- **Check parameters**: `./snaprepair --help`

## Example Output

```
SNAP-Guided Spectral Repair
============================
Input graph: ../as20graph.txt
Output prefix: demo_repair
Budget: 100 edges
Phi threshold: 0.200

Loading graph...
  Nodes: 6474, Edges: 13895

[...algorithm runs...]

Saving results...
  Edges saved to: demo_repair_edges.txt
  Modified graph saved to: demo_repair_graph.txt
  Final graph: Nodes=6474, Edges=13995

Total run time: 0.15s
```

---

**Need more help?** Check the full documentation in `ReadMe.txt` or run `./demo.sh` to see it in action!
