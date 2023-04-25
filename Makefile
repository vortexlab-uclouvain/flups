################################################################################
# @copyright Copyright (c) UCLouvain 2020
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
ARCH_FILE ?= make_arch/make.default

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
UNAME := $(shell uname)

#-----------------------------------------------------------------------------
CXX ?= mpic++
AR ?= ar

#-----------------------------------------------------------------------------
# Do not show GNU makefile command and info 
#.SILENT:
# git commit
#GIT_COMMIT ?= $(shell git describe --always --dirty)
GIT_COMMIT ?= $(shell git rev-parse --short HEAD)

PREFIX ?= ./
NAME := flups
# library naming
TARGET_LIB_ISR := build/lib$(NAME)_isr
TARGET_LIB_A2A := build/lib$(NAME)_a2a
TARGET_LIB_NB  := build/lib$(NAME)_nb

TARGET_LIB_DPREC_A2A := build/lib$(NAME)_dprec_a2a
TARGET_LIB_DPREC_NB  := build/lib$(NAME)_dprec_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build
KERNEL_DIR := ./kernel

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

#---- ACCFFT
#check if HAVE_ACCFFT
ifneq (,$(findstring -DHAVE_ACCFFT,$(CXXFLAGS)))
 ACCFFT_DIR ?= /usr/
 ACCFFT_INC ?= ${ACCFFT_DIR}/include
 ACCFFT_LIB ?= ${ACCFFT_DIR}/lib
 ACCFFT_LIBNAME ?= -laccfft
 INC += -I$(ACCFFT_INC)
 LIB += -L$(ACCFFT_LIB) $(ACCFFT_LIBNAME) -Wl,-rpath,$(ACCFFT_LIB)
endif

#-----------------------------------------------------------------------------
# LGF SPECIAL CASE
# by default the LGF kernel data is installed in the include directory
LGF_PATH=$(abspath $(PREFIX)/include)
DEF += -DKERNEL_PATH=${LGF_PATH}
LGF_DATA := $(wildcard $(KERNEL_DIR)/*.ker)

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)
API := $(wildcard $(SRC_DIR)/*.h)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ_ISR := $(SRC:%.cpp=$(OBJ_DIR)/isr_%.o)
OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)
OBJ_DPREC_A2A := $(SRC:%.cpp=$(OBJ_DIR)/dprec_a2a_%.o)
OBJ_DPREC_NB := $(SRC:%.cpp=$(OBJ_DIR)/dprec_nb_%.o)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)

################################################################################
# mandatory flags
M_FLAGS := -fPIC -DGIT_COMMIT=\"$(GIT_COMMIT)\" 

################################################################################
$(OBJ_DIR)/isr_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_ISR $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_NONBLOCK $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/dprec_nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_DPREC -DCOMM_NONBLOCK $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/dprec_a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_DPREC $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MMD -E $< -o $@

################################################################################
default: 
	@$(MAKE) info 
	@$(MAKE) lib_static

# for the validation, do a static lib
validation: install_static

# for the validation, do a static lib
test : install_static

# compile static and dynamic lib
all: lib_static lib_dynamic lib_darwin

all2all: $(TARGET_LIB_A2A).a $(TARGET_LIB_A2A).so

all2all_dprec: $(TARGET_LIB_DPREC_A2A).a $(TARGET_LIB_DPREC_A2A).so

nonblocking: $(TARGET_LIB_NB).a $(TARGET_LIB_NB).so $(TARGET_LIB_ISR).so $(TARGET_LIB_ISR).a

nonblocking_dprec: $(TARGET_LIB_DPREC_NB).a $(TARGET_LIB_DPREC_NB).so 

lib_static: $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a $(TARGET_LIB_ISR).a

lib_dynamic: $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so $(TARGET_LIB_ISR).so

ifeq ($(UNAME), Darwin)
lib_darwin: $(TARGET_LIB_A2A).dylib $(TARGET_LIB_NB).dylib $(TARGET_LIB_ISR).dylib
else
lib_darwin: 
endif

lib_static_deprec: $(TARGET_LIB_DPREC_A2A).a $(TARGET_LIB_DPREC_NB).a

lib_dynamic_deprec: $(TARGET_LIB_DPREC_A2A).so $(TARGET_LIB_DPREC_NB).so

ifeq ($(UNAME), Darwin)
lib_darwin_deprec: $(TARGET_LIB_DEPREC_A2A).dylib $(TARGET_LIB_DEPREC_NB).dylib
else
lib_darwin_deprec: 
endif

lib: lib_static

$(TARGET_LIB_ISR).so: $(OBJ_ISR)
	$(CXX) -shared $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).so: $(OBJ_A2A)
	$(CXX) -shared $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).so: $(OBJ_NB)
	$(CXX) -shared $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_DPREC_A2A).so: $(OBJ_DPREC_A2A)
	$(CXX) -shared $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_DPREC_NB).so: $(OBJ_DPREC_NB)
	$(CXX) -shared $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_ISR).a: $(OBJ_ISR)
	$(AR) rvs $(M_LFLAGS) $@  $^

$(TARGET_LIB_A2A).a: $(OBJ_A2A)
	$(AR) rvs $(M_LFLAGS) $@  $^

$(TARGET_LIB_NB).a: $(OBJ_NB)
	$(AR) rvs $(M_LFLAGS) $@  $^

$(TARGET_LIB_DPREC_A2A).a: $(OBJ_DPREC_A2A)
	$(AR) rvs $(M_LFLAGS) $@  $^

$(TARGET_LIB_DPREC_NB).a: $(OBJ_DPREC_NB)
	$(AR) rvs $(M_LFLAGS) $@  $^

# darwin specific target
darwin_dir := $(realpath $(PREFIX)/lib)
$(TARGET_LIB_ISR).dylib: $(OBJ_ISR)
	$(CXX) -dynamiclib -install_name $(darwin_dir)/$(notdir $@) $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).dylib: $(OBJ_A2A)
	$(CXX) -dynamiclib -install_name $(darwin_dir)/$(notdir $@)  $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).dylib: $(OBJ_NB)
	$(CXX) -dynamiclib -install_name $(darwin_dir)/$(notdir $@) $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_DPREC_A2A).dylib: $(OBJ_DPREC_A2A)
	$(CXX) -dynamiclib -install_name $(darwin_dir)/$(notdir $@) $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_DPREC_NB).dylib: $(OBJ_DPREC_NB)
	$(CXX) -dynamiclib -install_name $(darwin_dir)/$(notdir $@) $(LDFLAGS) $(M_LFLAGS) $^ -o $@ $(LIB)
# end of darwin targets

preproc: $(IN)

install_dynamic: lib_dynamic lib_darwin
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_ISR).so $(PREFIX)/lib
	@cp $(TARGET_LIB_A2A).so $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).so $(PREFIX)/lib
ifeq ($(UNAME), Darwin)
	@cp $(TARGET_LIB_ISR).dylib $(PREFIX)/lib
	@cp $(TARGET_LIB_A2A).dylib $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).dylib $(PREFIX)/lib
endif
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

install_static: lib_static 
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_ISR).a $(PREFIX)/lib
	@cp $(TARGET_LIB_A2A).a $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).a $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include


install_dprec_static: lib_static_deprec
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_DPREC_A2A).a $(PREFIX)/lib
	@cp $(TARGET_LIB_DPREC_NB).a $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

install_dprec_dynamic: lib_dynamic_deprec lib_darwin_deprec
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_DPREC_A2A).so $(PREFIX)/lib
	@cp $(TARGET_LIB_DPREC_NB).so $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include
# for a standard installation, do the dynamic link	
install: 
	@$(MAKE) info
	@$(MAKE) install_static
	@$(MAKE) install_dynamic

test:
	@echo $(SRC)

clean:
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TARGET_LIB_ISR).so $(TARGET_LIB_ISR).a
	@rm -f $(TARGET_LIB_A2A).so $(TARGET_LIB_A2A).a
	@rm -f $(TARGET_LIB_NB).so $(TARGET_LIB_NB).a
	@rm -f $(TARGET_LIB_DPREC_A2A).so $(TARGET_LIB_DPREC_A2A).a
	@rm -f $(TARGET_LIB_DPREC_NB).so $(TARGET_LIB_DPREC_NB).a

destroy:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.d
	@rm -f $(TARGET_LIB_ISR).so $(TARGET_LIB_ISR).a
	@rm -f $(TARGET_LIB_A2A).so $(TARGET_LIB_A2A).a
	@rm -f $(TARGET_LIB_NB).so $(TARGET_LIB_NB).a
	@rm -f $(TARGET_LIB_DPREC_A2A).so $(TARGET_LIB_DPREC_A2A).a
	@rm -f $(TARGET_LIB_DPREC_NB).so $(TARGET_LIB_DPREC_NB).a
	@rm -rf $(OBJ_DIR)/*
	@rm -rf include
	@rm -rf lib

#-------------------------------------------------------------------------------
# mentioning this target will export all the current variables to child-make processes
.EXPORT_ALL_VARIABLES:

.PHONY: info
info: 
	@$(MAKE) --file=info.mak
#-------------------------------------------------------------------------------

-include $(DEP)
