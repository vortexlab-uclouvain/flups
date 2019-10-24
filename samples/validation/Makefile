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
ARCH_FILE ?= make_arch/make.vagrant_intel
include ../../$(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
NAME := flups
# executable naming
TARGET_EXE_A2A := $(NAME)_validation_a2a
TARGET_EXE_NB := $(NAME)_validation_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FFTW
FFTW_INC ?= /usr/include
FFTW_LIB ?= /usr/lib
FFTW_LIBNAME ?= -lfftw3_omp -lfftw3
INC += -I$(FFTW_INC)
LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB)

#---- HDF5
HDF5_INC ?= /usr/include
HDF5_LIB ?= /usr/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB)

#---- FLUPS
FLUPS_INC ?= ../../include
FLUPS_LIB ?= ../../lib
INC += -I$(FLUPS_INC)

#---- METIS

INC+= -I$(METIS_INC)
LIB+= -L$(METIS_LIB) -lmetis  -Wl,-rpath,$(METIS_LIB)

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: all

all: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

all2all: $(TARGET_EXE_A2A)

nonblocking: $(TARGET_EXE_NB)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB) 

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB) 

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)

destroy:
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(OBJ_DIR)/*
	rm -rf include
	rm -rf lib

info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
