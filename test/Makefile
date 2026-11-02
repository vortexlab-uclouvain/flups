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
ARCH_FILE ?= make_arch/make.docker_clang_test
include ../$(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-------------------------------------------------------------------------------
# make predefines CXX (g++), so a plain ?= would never apply
ifeq ($(origin CXX),default)
CXX := mpic++
endif

#---- Detect nvcc (also behind an MPI wrapper) for the rpath syntax and the nvcc-only flags
# nvcc uses -Xlinker "-rpath=", while GNU/Clang use -Wl,-rpath,
IS_NVCC := $(findstring nvcc,$(shell $(CXX) --version 2>/dev/null))
ifneq (,$(IS_NVCC))
RPATH_FLAG = -Xlinker "-rpath=
RPATH_END = "
else
RPATH_FLAG = -Wl,-rpath,
RPATH_END =
endif

#-------------------------------------------------------------------------------
NAME := flups
# executable naming
TARGET_EXE_ISR := $(NAME)_test_isr
TARGET_EXE_A2A := $(NAME)_test_a2a
TARGET_EXE_NB := $(NAME)_test_nb

#-------------------------------------------------------------------------------
SRC_DIR := ./src
OBJ_DIR := ./build

## add the headers to the includes
INC := -I$(SRC_DIR)

#-------------------------------------------------------------------------------
# Required Dependencies
#---- FLUPS
FLUPS_DIR ?= ..
FLUPS_INC ?= $(FLUPS_DIR)/include
FLUPS_LIB ?= $(FLUPS_DIR)/lib
INC += -I$(FLUPS_INC)

#---- HDF5
HDF5_DIR ?= /usr
HDF5_INC ?= $(HDF5_DIR)/include
HDF5_LIB ?= $(HDF5_DIR)/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) $(RPATH_FLAG)$(HDF5_LIB)$(RPATH_END)

#---- GTEST
GTEST_DIR ?= /usr
GTEST_INC ?= $(GTEST_DIR)/include
GTEST_LIB ?= $(GTEST_DIR)/lib
GTEST_LIBNAME ?= -lgtest
INC += -I$(GTEST_INC)
LIB += -L$(GTEST_LIB) $(GTEST_LIBNAME) $(RPATH_FLAG)$(GTEST_LIB)$(RPATH_END)

#---- H3LPR
H3LPR_DIR ?= /usr
H3LPR_INC ?= $(H3LPR_DIR)/include
H3LPR_LIB ?= $(H3LPR_DIR)/lib
H3LPR_LIBNAME ?= -lh3lpr
INC += -I$(H3LPR_INC)
LIB += -L$(H3LPR_LIB) $(H3LPR_LIBNAME) $(RPATH_FLAG)$(H3LPR_LIB)$(RPATH_END)

#---- METIS (deprecated, see HAVE_METIS in the README)
ifneq (,$(findstring -DHAVE_METIS,$(CXXFLAGS) $(OPTS)))
METIS_DIR ?= /usr
METIS_INC ?= $(METIS_DIR)/include
METIS_LIB ?= $(METIS_DIR)/lib
INC += -I$(METIS_INC)
LIB += -L$(METIS_LIB) -lmetis $(RPATH_FLAG)$(METIS_LIB)$(RPATH_END)
endif

#-------------------------------------------------------------------------------
# Backend -- FFTW by default; the tests are host code, only the libraries are needed
ifeq (,$(findstring -DBACKEND_,$(OPTS)))
OPTS += -DBACKEND_FFTW
endif

#---- FFTW
ifneq (,$(findstring -DBACKEND_FFTW,$(OPTS)))
FFTW_DIR ?= /usr
FFTW_INC ?= $(FFTW_DIR)/include
FFTW_LIB ?= $(FFTW_DIR)/lib
FFTW_LIBNAME ?= -lfftw3
INC += -I$(FFTW_INC)
LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) $(RPATH_FLAG)$(FFTW_LIB)$(RPATH_END)
endif

#---- CUDA
# the stubs are only for linking on a node without driver, never in the rpath
ifneq (,$(findstring -DBACKEND_CUDA,$(OPTS)))
CUDA_DIR ?= /usr/local/cuda
CUDA_INC ?= $(CUDA_DIR)/include
CUDA_LIB ?= $(CUDA_DIR)/lib64
CUDA_LIBNAME ?= -lcudart -lcufft -lcublas
INC += -I$(CUDA_INC)
LIB += -L$(CUDA_LIB) $(CUDA_LIBNAME) $(RPATH_FLAG)$(CUDA_LIB)$(RPATH_END) -L$(CUDA_LIB)/stubs
ifneq (,$(IS_NVCC))
CXXFLAGS += --forward-unknown-to-host-compiler
LDFLAGS += --forward-unknown-to-host-compiler
endif
endif

#---- HIP
ifneq (,$(findstring -DBACKEND_HIP,$(OPTS)))
HIP_DIR ?= /opt/rocm
HIP_INC ?= $(HIP_DIR)/include
HIP_LIB ?= $(HIP_DIR)/lib
HIP_LIBNAME ?= -lrocblas -lhipfft -lamdhip64
INC += -I$(HIP_INC)
LIB += -L$(HIP_LIB) $(HIP_LIBNAME) $(RPATH_FLAG)$(HIP_LIB)$(RPATH_END)
endif

#-------------------------------------------------------------------------------
## sources and headers
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
.PHONY: default all all2all nonblocking
default: all

all: info $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_ISR)

all2all: $(TARGET_EXE_A2A)

nonblocking: $(TARGET_EXE_NB) $(TARGET_EXE_ISR)

#-------------------------------------------------------------------------------
$(TARGET_EXE_ISR): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -l$(NAME)_isr $(RPATH_FLAG)$(FLUPS_LIB)$(RPATH_END) $(LIB)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -l$(NAME)_a2a $(RPATH_FLAG)$(FLUPS_LIB)$(RPATH_END) $(LIB)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -l$(NAME)_nb $(RPATH_FLAG)$(FLUPS_LIB)$(RPATH_END) $(LIB)

#-------------------------------------------------------------------------------
.PHONY: clean destroy
clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_ISR) $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

destroy: clean
	rm -f $(OBJ_DIR)/*.d

#-------------------------------------------------------------------------------
.PHONY: info
info:
	$(info compiler = $(shell $(CXX) --version))
	$(info compil. flags = $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD)
	$(info linker flags = $(LDFLAGS))
	$(info libraries = -L$(FLUPS_LIB) $(RPATH_FLAG)$(FLUPS_LIB)$(RPATH_END) $(LIB))
	$(info using arch file = $(ARCH_FILE))
	$(info ------------)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info ------------)

-include $(DEP)

# mkdir the needed dirs
$(shell mkdir -p $(OBJ_DIR))
