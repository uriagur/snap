#!/bin/bash
# Compare NCP plots of original and repaired graphs

if [ $# -lt 2 ]; then
    echo "Usage: $0 <original_graph> <repaired_graph_prefix>"
    echo "Example: $0 ../as20graph.txt test_repair"
    exit 1
fi

ORIGINAL_GRAPH=$1
REPAIRED_PREFIX=$2
REPAIRED_GRAPH="${REPAIRED_PREFIX}_graph.txt"

# Check if files exist
if [ ! -f "$ORIGINAL_GRAPH" ]; then
    echo "Error: Original graph file not found: $ORIGINAL_GRAPH"
    exit 1
fi

if [ ! -f "$REPAIRED_GRAPH" ]; then
    echo "Error: Repaired graph file not found: $REPAIRED_GRAPH"
    echo "Did you run snaprepair first?"
    exit 1
fi

# Navigate to ncpplot directory
cd ../ncpplot

echo "Generating NCP plot for ORIGINAL graph..."
./ncpplot -i:"../../$ORIGINAL_GRAPH" -o:original -d:"Original Graph" -kmin:10 -kmax:100000 -c:10

echo ""
echo "Generating NCP plot for REPAIRED graph..."
./ncpplot -i:../snaprepair/"$REPAIRED_GRAPH" -o:repaired -d:"SNAP-Repaired Graph" -kmin:10 -kmax:100000 -c:10

echo ""
echo "Done! NCP plots generated:"
echo "  Original: ncpplot/ncp.original.*.png"
echo "  Repaired: ncpplot/ncp.repaired.*.png"
echo ""
echo "You can compare these plots to see the effect of spectral repair."
