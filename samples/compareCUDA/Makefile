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
ARCH_FILE ?= make_arch/make.docker_clang
include ../../$(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
NAME := flups_gpu_vs_the_world
# executable naming
TARGET_EXE_ISR := $(NAME)_isr
TARGET_EXE_A2A := $(NAME)_a2a
TARGET_EXE_NB := $(NAME)_nb
TARGET_EXE_DPREC_A2A := $(NAME)_validation_dprec_a2a
TARGET_EXE_DPREC_NB := $(NAME)_validation_dprec_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FLUPS
FLUPS_INC ?= ../../include
FLUPS_LIB ?= ../../lib
INC += -I$(FLUPS_INC)

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

#---- H3LPR
H3LPR_INC ?= /usr/include
H3LPR_LIB ?= /usr/lib
H3LPR_LIBNAME ?= -lh3lpr
INC += -I$(H3LPR_INC)
LIB += -L$(H3LPR_LIB) $(H3LPR_LIBNAME) -Wl,-rpath,$(H3LPR_LIB)

#---- HeFFTE
ifneq (,$(findstring -DHAVE_HEFFTE,$(CXXFLAGS)))
	HEFFTE_INC ?= /usr/include
	HEFFTE_LIB ?= /usr/lib
	HEFFTE_LIBNAME ?= -lheffte
	INC += -I$(HEFFTE_INC)
	LIB += -L$(HEFFTE_LIB) $(HEFFTE_LIBNAME) -Wl,-rpath,$(HEFFTE_LIB)
endif

#-----------------------------------------------------------------------------
# Pick a backend -- by default, we pick FFTW
ifeq (,$(findstring -DBACKEND_,$(OPTS)))
	OPTS += -DBACKEND_FFTW
endif

#---- Use FFTW
ifneq (,$(findstring -DBACKEND_FFTW,$(OPTS)))
    FFTW_INC ?= /usr/include
    FFTW_LIB ?= /usr/lib
    FFTW_LIBNAME ?= -lfftw3_omp -lfftw3
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
	ifneq (,$(findstring nvcc,$(shell $(CXX) --version)))
		CXXFLAGS += --forward-unknown-to-host-compiler
		LDFLAGS  += --forward-unknown-to-host-compiler
	endif
endif

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(HEAD)
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_ISR)

all: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_ISR)


all2all: $(TARGET_EXE_A2A)

nonblocking: $(TARGET_EXE_NB) $(TARGET_EXE_ISR)

nonblockingisr: $(TARGET_EXE_ISR)

deprec: $(TARGET_EXE_DPREC_A2A) $(TARGET_EXE_DPREC_NB)

nonblocking_deprec: $(TARGET_EXE_DPREC_NB)

all2all_deprec: $(TARGET_EXE_DPREC_A2A)

$(TARGET_EXE_ISR): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_isr -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_DPREC_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_dprec_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_DPREC_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_dprec_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_ISR)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(TARGET_EXE_DPREC_A2A)
	rm -f $(TARGET_EXE_DPREC_NB)

destroy: clean
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
