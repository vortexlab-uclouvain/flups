################################################################################
# @copyright Copyright © UCLouvain 2020
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
ARCH_FILE ?= make_arch/make.osx_gcc

include ../../$(ARCH_FILE)

FLUPSDIR := ../../

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := solve_advanced
# executable naming
TARGET_EXE_A2A := $(NAME)_a2a
TARGET_EXE_NB := $(NAME)_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

#-----------------------------------------------------------------------------
#---- FLUPS
INC += -I$(FLUPSDIR)/include
LIB_A2A = -L$(FLUPSDIR)/lib -lflups_a2a
LIB_NB  = -L$(FLUPSDIR)/lib -lflups_nb 

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

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.c))
# HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
# DEP := $(SRC:%.c=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.c=$(OBJ_DIR)/%.o)

################################################################################

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CCFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

all: default

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB_A2A) $(LIB)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB_NB) $(LIB)

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
	
destroy:
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(OBJ_DIR)/*
	
info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	# $(info DEP = $(DEP))

# -include $(DEP)
