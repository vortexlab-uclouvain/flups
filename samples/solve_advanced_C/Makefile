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

#---- HDF5
HDF5_INC ?= /usr/include
HDF5_LIB ?= /usr/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB)

#---- H3LPR
H3LPR_INC ?= /usr/include
H3LPR_LIB ?= /usr/lib
H3LPR_LIBNAME ?= -lh3lpr
INC += -I$(H3LPR_INC)
LIB += -L$(H3LPR_LIB) $(H3LPR_LIBNAME) -Wl,-rpath,$(H3LPR_LIB)

#-----------------------------------------------------------------------------
# Pick a backend -- by default, we pick FFTW
ifeq (,$(findstring -DBACKEND_,$(OPTS)))
	OPTS += -DBACKEND_FFTW
endif

#---- Use FFTW
ifneq (,$(findstring -DBACKEND_FFTW,$(OPTS)))
    FFTW_INC ?= /usr/include
    FFTW_LIB ?= /usr/lib
    FFTW_LIBNAME ?= -lfftw3
    INC += -I$(FFTW_INC)
    LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB)
endif

#---- Use CUDA
ifneq (,$(findstring -DBACKEND_CUDA,$(OPTS)))
	CUDA_INC ?= /usr/local/cuda/include
	CUDA_LIB ?= /usr/local/cuda/lib64
	CUDA_LIBNAME ?= -lcudart -lcufft -lcublas
	INC += -I$(CUDA_INC)
	LIB += -L$(CUDA_LIB) $(CUDA_LIBNAME) -Wl,-rpath,$(CUDA_LIB)
	LIB += -L$(CUDA_LIB)/stubs -Wl,-rpath,$(CUDA_LIB)/stubs
endif

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
