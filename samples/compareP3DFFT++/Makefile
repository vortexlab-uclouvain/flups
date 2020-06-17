################################################################################
# @copyright Copyright © UCLouvain 2020
# 
# FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
# 
# Copyright (C) <2020> <Universite catholique de Louvain (UCLouvain), Belgique>
# 
# List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#  http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.osx_gcc
#DYLD_LIBRARY_PATH=/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver/lib OMPI_MCA_btl=self,tcp OMPI_MCA_rmaps_base_oversubscribe=yes OMP_NUM_THREADS=1 mpirun -np 1 --tag-output --oversubscribe ./flups_vs_p3dfft 

include ../../$(ARCH_FILE)
#include $(ARCH_FILE)

FLUPSDIR := ../../
#/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver
P3DFFTDIR++ := ./

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := flups_vs_p3dfft++
# executable naming
TARGET_EXE_A2A := $(NAME)_a2a
TARGET_EXE_NB := $(NAME)_nb
TARGET_EXE_A2A_SKIP := $(NAME)_a2a_noP3D
TARGET_EXE_NB_SKIP := $(NAME)_nb_noP3D
#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FLUPS
INC += -I$(FLUPSDIR)/include -I$(FLUPSDIR)/src 
LIB_A2A = -L$(FLUPSDIR)/lib -lflups_a2a
LIB_NB  = -L$(FLUPSDIR)/lib -lflups_nb 

#---- P3DFFT++
INC += -I$(P3DFFTDIR++)/include -DHAVE_CONFIG_H
LIB += -L$(P3DFFTDIR++)/lib -lp3dfft.3

#---- FFTW
FFTW_INC ?= /usr/include
FFTW_LIB ?= /usr/lib
FFTW_LIBNAME ?= -lfftw3_omp -lfftw3 -lfftw3f 
INC += -I$(FFTW_INC)
LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB)

#---- HDF5
HDF5_INC ?= /usr/include
HDF5_LIB ?= /usr/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB)

#---- METIS
#check if HAVE_METIS
ifneq (,$(findstring -DHAVE_METIS,$(CXXFLAGS)))
	METIS_INC ?= /usr/include
	METIS_LIB ?= /usr/lib
	INC+= -I$(METIS_INC)
	LIB+= -L$(METIS_LIB) -lmetis  -Wl,-rpath,$(METIS_LIB)
endif

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
OBJ_SKIP := $(SRC:%.cpp=$(OBJ_DIR)/noP3D_%.o)

################################################################################

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/noP3D_%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -DSKIP_P3D -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

skip: $(TARGET_EXE_A2A_SKIP) $(TARGET_EXE_NB_SKIP)

all: default

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_A2A)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_NB)

$(TARGET_EXE_A2A_SKIP): $(OBJ_SKIP)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_A2A)

$(TARGET_EXE_NB_SKIP): $(OBJ_SKIP)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_NB)

install: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)
	mkdir -p $(PREFIX)/bin
	cp $(TARGET_EXE_A2A) $(PREFIX)/bin
	cp $(TARGET_EXE_NB) $(PREFIX)/bin

test:
	@echo $(SRC)
	@echo $(INC)
	@echo $(LIB)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(TARGET_EXE_A2A_SKIP)
	rm -f $(TARGET_EXE_NB_SKIP)
	
destroy: clean
	rm -f $(OBJ_DIR)/*
	
info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info OBJ_SKIP = $(OBJ_SKIP))
	$(info DEP = $(DEP))

-include $(DEP)
