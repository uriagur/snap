#!/bin/bash
#
# Build and run the Dynamic Repair algorithm
#

set -e  # Exit on error

echo "=================================="
echo "Dynamic Repair Algorithm - Builder"
echo "=================================="
echo ""

# Check if we're in the right directory
if [ ! -f "test_dynamic_repair.cpp" ]; then
    echo "Error: Must run from examples/snaprepair directory"
    exit 1
fi

# Build Snap.o if needed
echo "Step 1: Building SNAP core library..."
if [ ! -f "../../snap-core/Snap.o" ]; then
    echo "  Building snap-core..."
    make -C ../../snap-core
else
    echo "  snap-core already built (Snap.o exists)"
fi

# Build test_dynamic_repair
echo ""
echo "Step 2: Building test_dynamic_repair..."
g++ -std=c++98 -Wall -O3 -DNDEBUG -fopenmp \
    -I../../glib-core -I../../snap-core -I../../snap-adv \
    test_dynamic_repair.cpp RepairStrategies.cpp DynamicRepair.cpp \
    ../../snap-adv/ncp.cpp -o test_dynamic_repair \
    ../../snap-core/Snap.o -lrt

if [ -f "test_dynamic_repair" ]; then
    echo "  ✓ Build successful: test_dynamic_repair"
else
    echo "  ✗ Build failed"
    exit 1
fi

echo ""
echo "=================================="
echo "Build Complete!"
echo "=================================="
echo ""
echo "To run the algorithm:"
echo "  ./test_dynamic_repair -i:../as20graph.txt -b:100 -t:0.2"
echo ""
echo "Quick test (small budget):"
echo "  ./test_dynamic_repair -i:../as20graph.txt -b:10 -t:0.2 -kmin:10 -kmax:50 -c:2"
echo ""
echo "For help and parameters:"
echo "  ./test_dynamic_repair"
echo ""
echo "For detailed documentation:"
echo "  cat QUICKSTART_DYNAMIC_REPAIR.md"
echo ""
