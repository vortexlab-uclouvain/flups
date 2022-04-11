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
ARCH_FILE ?= make_arch/make.default

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
# Do not show GNU makefile command and info 
#.SILENT:

PREFIX ?= ./
NAME := flups
# library naming
TARGET_LIB_A2A := build/lib$(NAME)_a2a
TARGET_LIB_NB  := build/lib$(NAME)_nb
TARGET_LIB_AGGRESSIVE  := build/lib$(NAME)_agg

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
OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)
OBJ_AGGRESSIVE := $(SRC:%.cpp=$(OBJ_DIR)/mpiaggressive_%.o)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)

################################################################################
$(OBJ_DIR)/nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_NONBLOCK $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/mpiaggressive_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DMPI_AGGRESSIVE $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD -E $< -o $@

################################################################################
default: 
	@$(MAKE) info 
	@$(MAKE) lib_static 

# for the validation, do a static lib
validation: install_static

# for the validation, do a static lib
test : install_static

# compile static and dynamic lib
all: lib_static lib_dynamic lib_agressive

all2all: $(TARGET_LIB_A2A).a $(TARGET_LIB_A2A).so

nonblocking: $(TARGET_LIB_NB).a $(TARGET_LIB_NB).so

aggressive: $(TARGET_LIB_AGGRESSIVE).a $(TARGET_LIB_AGGRESSIVE).so

lib_static: $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a $(TARGET_LIB_AGGRESSIVE).a

lib_dynamic: $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so $(TARGET_LIB_AGGRESSIVE).so

lib: lib_static

$(TARGET_LIB_A2A).so: $(OBJ_A2A)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).so: $(OBJ_NB)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_AGGRESSIVE).so: $(OBJ_AGGRESSIVE)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).a: $(OBJ_A2A)
	ar rvs $@  $^  

$(TARGET_LIB_NB).a: $(OBJ_NB)
	ar rvs $@  $^  

$(TARGET_LIB_AGGRESSIVE).a: $(OBJ_AGGRESSIVE)
	ar rvs $@  $^  

preproc: $(IN)

install_dynamic: lib_dynamic
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_A2A).so $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).so $(PREFIX)/lib
	@cp $(TARGET_LIB_AGGRESSIVE).so $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

install_static: lib_static 
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_A2A).a $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).a $(PREFIX)/lib
	@cp $(TARGET_LIB_AGGRESSIVE).a $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

# for a standard installation, do the dynamic link	
install: 
	@$(MAKE) info
	@$(MAKE) install_static

test:
	@echo $(SRC)

clean:
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TARGET_LIB_A2A)*
	@rm -f $(TARGET_LIB_NB)*
	@rm -f $(TARGET_LIB_AGGRESSIVE)*

destroy:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.d
	@rm -rf $(TARGET_LIB_A2A)
	@rm -rf $(TARGET_LIB_NB)
	@rm -rf $(TARGET_LIB_AGGRESSIVE)
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
