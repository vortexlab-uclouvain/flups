################################################################################
# @copyright Copyright © UCLouvain 2019
# 
# FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
# 
# Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
# 
# List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
# 
# This program (FLUPS) is free software: 
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program (see COPYING file).  If not, 
# see <http://www.gnu.org/licenses/>.
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= ../make_arch/make.osx_gcc
#DYLD_LIBRARY_PATH=/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver/lib OMPI_MCA_btl=self,tcp OMPI_MCA_rmaps_base_oversubscribe=yes OMP_NUM_THREADS=1 mpirun -np 1 --tag-output --oversubscribe ./flups_vs_p3dfft 

include $(ARCH_FILE)

#CXXFLAGS := -g -Wall -O0 --debug -std=c++17 -DVERBOSE 

FLUPSDIR := ../
#/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver
P3DFFTDIR := ./

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := flups_vs_p3dfft
# executable naming
TARGET_EXE := $(NAME)
# TARGET_EXE_A2A := $(NAME)_validation_a2a
# TARGET_EXE_NB := $(NAME)_validation_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FLUPS
INC += -I$(FLUPSDIR)/include
LIB += -L$(FLUPSDIR)/lib -lflups_a2a -Wl,-rpath,$(FLUPSDIR)/lib

#---- P3DFFT
INC += -I$(P3DFFTDIR)/include
LIB += -L$(P3DFFTDIR)/lib -lp3dfft -Wl,-rpath,$(P3DFFTDIR)/lib -lmpi_mpifh -lgfortran 

#---- FFTW
INC += -I$(FFTWDIR)/include
LIB += -L$(FFTWDIR)/lib -lfftw3  -Wl,-rpath,$(FFTWDIR)/lib

#---- HDF5
HDF5LIB ?= -L$(HDF5DIR)/lib -lhdf5 -Wl,-rpath,$(HDF5DIR)/lib
HDF5INC ?= -I$(HDF5DIR)/include
INC += $(HDF5INC)
LIB += $(HDF5LIB)

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
# OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
# OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)

################################################################################

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE)

all: default

$(TARGET_EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

install: $(TARGET_EXE)
	mkdir -p $(PREFIX)/bin
	cp $(TARGET_EXE) $(PREFIX)/bin

test:
	@echo $(SRC)
	@echo $(INC)
	@echo $(LIB)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE)
	
destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET_EXE)
	rm -f $(OBJ_DIR)/*
	
info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
