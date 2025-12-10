#
#	configuration variables for the example

## Main application file
MAIN = ncpplot
DEPH = $(EXSNAPADV)/ncp.h
DEPCPP = $(EXSNAPADV)/ncp.cpp

## Debug flags
# Uncomment the next line to enable debugging
# DEBUG = 1

ifdef DEBUG
    CXXFLAGS := $(filter-out -O3 -DNDEBUG,$(CXXFLAGS))
    CXXFLAGS += -g -ggdb -O0 -DDEBUG
endif

