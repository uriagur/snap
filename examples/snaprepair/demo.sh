#!/bin/bash
# Demo script for SNAP-Guided Spectral Repair
# This script demonstrates the complete workflow

set -e  # Exit on error

echo "================================="
echo "SNAP-Guided Spectral Repair Demo"
echo "================================="
echo ""

# Configuration
GRAPH_FILE="../musae_facebook_edges.txt"
OUTPUT_PREFIX="musae_facebook_edges_repair10000"
BUDGET=10000
PHI_THRESHOLD=0.12
KMIN=10
KMAX=100000

echo "Configuration:"
echo "  Input Graph: $GRAPH_FILE"
echo "  Output Prefix: $OUTPUT_PREFIX"
echo "  Edge Budget: $BUDGET"
echo "  Phi Threshold: $PHI_THRESHOLD"
echo "  K Range: [$KMIN, $KMAX]"
echo ""

# Step 1: Run the repair algorithm
echo "Step 1: Running SNAP-Guided Spectral Repair..."
echo "----------------------------------------------"
./snaprepair \
    -i:"$GRAPH_FILE" \
    -o:"$OUTPUT_PREFIX" \
    -b:$BUDGET \
    -t:$PHI_THRESHOLD \
    -kmin:$KMIN \
    -kmax:$KMAX \
    -c:10 \
    -ncp:F

echo ""
echo "Step 2: Analyzing Results..."
echo "-----------------------------"

# Count original and modified edges
ORIG_EDGES=$(wc -l < "$GRAPH_FILE")
MOD_EDGES=$(wc -l < "${OUTPUT_PREFIX}_graph.txt")
ADDED_EDGES=$(tail -n +4 "${OUTPUT_PREFIX}_edges.txt" | wc -l)

echo "Original graph edges: $ORIG_EDGES"
echo "Modified graph edges: $MOD_EDGES"
echo "Edges added: $ADDED_EDGES"
echo ""

echo "Sample of added edges:"
tail -n +4 "${OUTPUT_PREFIX}_edges.txt" | head -10
echo "... (see ${OUTPUT_PREFIX}_edges.txt for full list)"
echo ""

echo "Step 3: Generating NCP Plots for Comparison..."
echo "------------------------------------------------"

# Build ncpplot if needed
if [ ! -f ../ncpplot/ncpplot ]; then
    echo "Building ncpplot..."
    cd ../ncpplot && make && cd ../snaprepair
fi

# Generate NCP plot for repaired graph
echo "Generating NCP for repaired graph..."
cd ../ncpplot
./ncpplot -i:../snaprepair/"${OUTPUT_PREFIX}_graph.txt" -k:F
cd ../snaprepair

echo ""
echo "================================="
echo "Demo Complete!"
echo "================================="
echo ""
echo "Generated Files:"
echo "  1. ${OUTPUT_PREFIX}_edges.txt - List of added edges"
echo "  2. ${OUTPUT_PREFIX}_graph.txt - Modified graph"
echo "  3. ../ncpplot/ncp.demo_original.*.png - Original NCP plot"
echo "  4. ../ncpplot/ncp.demo_repaired.*.png - Repaired NCP plot"
echo ""
echo "What to look for:"
echo "  - The NCP curve should be flatter (higher at small k values)"
echo "  - 'Dips' or 'valleys' in the original should be reduced"
echo "  - Overall conductance should be more uniform across cluster sizes"
echo ""
echo "To view the plots, open the PNG files in ../ncpplot/"
