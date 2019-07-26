################################################################################
# 
#-----------------------------------------------------------------------------
# machine dependent file: define INC_DIR ; LIB_DIR ; LIB
#include make_arch/make.VVPM
#include make.zenobe

TARGET := poisson

#-----------------------------------------------------------------------------
# COMPILER AND OPT/DEBUG FLAGS
CC = mpiicc
CXX = mpiicpc

CXXFLAGS := -O3 -g -DNDEBUG -stdc++11 -qopt-report=5
# CXXFLAGS := -g -Wall -O0 -traceback -ftrapuv -debug all -stdc++11 -DVERBOSE 
# CXXFLAGS := -g -Wall -O0 -traceback -ftrapuv -debug all -stdc++11

# define some options
DEF :=

# linker specific flags
LDFLAGS := 

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------

#---- FFTW
FFTWDIR  := /vagrant/soft/fftw-3.3.8-intel_2019.04
INC += -I$(FFTWDIR)/include
LIB += -L$(FFTWDIR)/lib -lfftw3_mpi -lfftw3

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################

default: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIB)

test:
	@echo $(SRC)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET)

destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET)
	rm -f $(OBJ_DIR)/*

info:
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
