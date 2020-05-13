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
ARCH_FILE ?= make_arch/make.generic

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := flups
# library naming
TARGET_LIB_A2A := build/lib$(NAME)_a2a
TARGET_LIB_NB  := build/lib$(NAME)_nb

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
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)

################################################################################
$(OBJ_DIR)/nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) -DCOMM_NONBLOCK $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -E $< -o $@

################################################################################
default: lib_static 

# for the validation, do a static lib
validation: install_static

# compile static and dynamic lib
all: lib_static lib_dynamic

all2all: $(TARGET_LIB_A2A).a $(TARGET_LIB_A2A).so

nonblocking: $(TARGET_LIB_NB).a $(TARGET_LIB_NB).so

lib_static: info $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a

lib_dynamic: info $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so

lib: lib_static

$(TARGET_LIB_A2A).so: $(OBJ_A2A)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).so: $(OBJ_NB)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).a: $(OBJ_A2A)
	ar rvs $@  $^  

$(TARGET_LIB_NB).a: $(OBJ_NB)
	ar rvs $@  $^  

preproc: $(IN)

install_dynamic: lib_dynamic
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_A2A).so $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).so $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

install_static: lib_static 
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@cp $(TARGET_LIB_A2A).a $(PREFIX)/lib
	@cp $(TARGET_LIB_NB).a $(PREFIX)/lib
	@cp $(API) $(PREFIX)/include
	@cp $(LGF_DATA) $(PREFIX)/include

# for a standard installation, do the dynamic link	
install: info install_static

test:
	@echo $(SRC)

clean:
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TARGET_LIB_A2A)*
	@rm -f $(TARGET_LIB_NB)*

destroy:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.d
	@rm -rf $(TARGET_LIB_A2A)
	@rm -rf $(TARGET_LIB_NB)
	@rm -rf $(OBJ_DIR)/*
	@rm -rf include
	@rm -rf lib

info: logo
	$(info prefix = $(PREFIX)/lib )
	$(info compiler = $(shell $(CXX) --version))
	$(info compil. flags = $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD)
	$(info linker flags = -shared $(LDFLAGS))
	$(info using arch file = $(ARCH_FILE) )
	$(info LGF path = $(LGF_PATH) )
	$(info ------------)
	$(info FFTW:)
	$(info - include: -I$(FFTW_INC) )
	$(info - lib: -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB))
	$(info ------------)
	$(info HDF5:)
	$(info - include: -I$(HDF5_INC) )
	$(info - lib: -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB))
	$(info ------------)
	$(info LIST OF OBJECTS:)
	$(info - SRC = $(SRC))
	$(info - OBJ A2A = $(OBJ_A2A))
	$(info - OBJ NB = $(OBJ_NB))
	$(info - DEP = $(DEP))
	$(info - LGF_DATA = $(LGF_DATA))
	$(info ------------)

.NOTPARALLEL: logo

logo: 
	@echo "----------------------------------------------------"
	@echo "    ______   _        _    _   _____     _____       "
	@echo "   |  ____| | |      | |  | | |  __ \   / ____|     "
	@echo "   | |__    | |      | |  | | | |__) | | (___       "
	@echo "   |  __|   | |      | |  | | |  ___/   \___ \      "
	@echo "   | |      | |____  | |__| | | |       ____) |     "
	@echo "   |_|      |______|  \____/  |_|      |_____/      "
	@echo "                                                    "
	@echo "                                                    "
	@echo "    	(C) UCLouvain - GPLv3                          "
	@echo "----------------------------------------------------"

-include $(DEP)
