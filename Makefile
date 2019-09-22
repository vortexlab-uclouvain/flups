################################################################################
# @copyright Copyright © UCLouvain 2019
# 
# FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
# 
# Copyright (C) <2019> <Université catholique de Louvain (UCLouvain), Belgique>
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
include make_arch/make.vagrant_intel

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
NAME := flups
TARGET_EXE := $(NAME)
TARGET_LIB := build/lib$(NAME).so

PREFIX ?= ./

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FFTW
INC += -I$(FFTWDIR)/include
LIB += -L$(FFTWDIR)/lib -lfftw3_omp -lfftw3  -Wl,-rpath,$(FFTWDIR)/lib

#---- HDF5
INC += -I$(HDF5DIR)/include
LIB += -L$(HDF5DIR)/lib -lhdf5 -Wl,-rpath,$(HDF5DIR)/lib

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE)

$(TARGET_EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB): $(OBJ)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

install: $(TARGET_LIB)
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $(TARGET_LIB) $(PREFIX)/lib
	cp $(HEAD) $(PREFIX)/include

test:
	@echo $(SRC)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_LIB)

destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_LIB)
	rm -f $(OBJ_DIR)/*
	rm -f include/*
	rm -f lib/*

info:
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
