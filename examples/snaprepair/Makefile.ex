#
#	configuration variables for the example

## Main application file
MAIN = snaprepair
DEPH = $(EXSNAPADV)/ncp.h ClusterManager.hpp RepairStrategies.h DynamicRepair.h
DEPCPP = $(EXSNAPADV)/ncp.cpp RepairStrategies.cpp DynamicRepair.cpp

