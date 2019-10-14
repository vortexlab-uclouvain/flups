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
ARCH_FILE ?= make_arch/make.osx_gcc

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := flups
# executable naming
TARGET_EXE := $(NAME)_validation
TARGET_EXE_A2A := $(NAME)_validation_a2a
TARGET_EXE_NB := $(NAME)_validation_nb
# library naming
TARGET_LIB_A2A := build/lib$(NAME)_a2a
TARGET_LIB_NB  := build/lib$(NAME)_nb

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
HDF5LIB ?= -L$(HDF5DIR)/lib -lhdf5 -Wl,-rpath,$(HDF5DIR)/lib
HDF5INC ?= -I$(HDF5DIR)/include
INC += $(HDF5INC)
LIB += $(HDF5LIB)

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)
API := $(wildcard $(SRC_DIR)/*.h)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)

################################################################################
$(OBJ_DIR)/nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) -DCOMM_NONBLOCK $(INC) $(DEF) -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE)

all: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) lib_static

validation: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

all2all: $(TARGET_EXE_A2A)

nonblocking: $(TARGET_EXE_NB)

lib_static: $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a

lib_dynamic: $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so

lib: lib_static

$(TARGET_EXE): $(OBJ_A2A)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_EXE_A2A): $(OBJ_A2A)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_EXE_NB): $(OBJ_NB)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).so: $(OBJ_A2A)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).so: $(OBJ_NB)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).a: $(OBJ_A2A)
	ar rvs $@  $^  

$(TARGET_LIB_NB).a: $(OBJ_NB)
	ar rvs $@  $^  

install_dynamic: lib_dynamic
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $(TARGET_LIB_A2A).so $(PREFIX)/lib
	cp $(TARGET_LIB_NB).so $(PREFIX)/lib
	cp $(HEAD) $(PREFIX)/include

install_static: lib_static
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $(TARGET_LIB_A2A).a $(PREFIX)/lib
	cp $(TARGET_LIB_NB).a $(PREFIX)/lib
	cp $(HEAD) $(PREFIX)/include
	cp $(SRC_DIR)/flups.h $(PREFIX)/include/flups.h
	
install: install_static

test:
	@echo $(SRC)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(TARGET_LIB_A2A)*
	rm -f $(TARGET_LIB_NB)*

destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(TARGET_LIB_A2A)*
	rm -f $(TARGET_LIB_NB)*
	rm -f $(OBJ_DIR)/*
	rm -rf include
	rm -rf lib

info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ_A2A))
	$(info OBJ = $(OBJ_NB))
	$(info DEP = $(DEP))

-include $(DEP)
