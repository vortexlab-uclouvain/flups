################################################################################
# @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
#      See LICENSE file in top-level directory
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.osx_gcc
#DYLD_LIBRARY_PATH=/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver/lib OMPI_MCA_btl=self,tcp OMPI_MCA_rmaps_base_oversubscribe=yes OMP_NUM_THREADS=1 mpirun -np 1 --tag-output --oversubscribe ./flups_vs_p3dfft 

include ../../$(ARCH_FILE)
#include $(ARCH_FILE)

#CXXFLAGS := -g -Wall -O0 --debug -std=c++17 -DVERBOSE 

FLUPSDIR := ../../
#/Users/DeeGee/Documents/PhD/devel/mix-poisson-solver
P3DFFTDIR := ./

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := flups_vs_accfft
# executable naming
TARGET_EXE_ISR := $(NAME)_isr
TARGET_EXE_A2A := $(NAME)_a2a
TARGET_EXE_NB := $(NAME)_nb
TARGET_EXE_DPREC_A2A := $(NAME)_dprec_a2a
TARGET_EXE_DPREC_NB := $(NAME)_dprec_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./
OBJ_DIR := ./build

#-----------------------------------------------------------------------------
#---- FLUPS
INC += -I$(FLUPSDIR)/include -I$(FLUPSDIR)/src 
LIB_ISR = -L$(FLUPSDIR)/lib -lflups_isr
LIB_A2A = -L$(FLUPSDIR)/lib -lflups_a2a
LIB_NB  = -L$(FLUPSDIR)/lib -lflups_nb 
LIB_DPREC_A2A = -L$(FLUPSDIR)/lib -lflups_dprec_a2a
LIB_DPREC_NB  = -L$(FLUPSDIR)/lib -lflups_dprec_nb 

#---- ACCFFT
ACCFFT_DIR ?= /usr/
ACCFFT_INC ?= ${ACCFFT_DIR}/include
ACCFFT_LIB ?= ${ACCFFT_DIR}/lib
ACCFFT_LIBNAME ?= -laccfft
INC += -I$(ACCFFT_INC)
LIB += -L$(ACCFFT_LIB) $(ACCFFT_LIBNAME) -Wl,-rpath,$(ACCFFT_LIB)

#---- FFTW
FFTW_INC ?= /usr/include
FFTW_LIB ?= /usr/lib
FFTW_LIBNAME ?= -lfftw3_omp -lfftw3 -lfftw3_threads
INC += -I$(FFTW_INC)
LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB)

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

#---- ACCFFT
ACCFFT_DIR ?= /usr/
ACCFFT_INC ?= ${ACCFFT_DIR}/include
ACCFFT_LIB ?= ${ACCFFT_DIR}/lib
ACCFFT_LIBNAME ?= -laccfft
INC += -I$(ACCFFT_INC)
LIB += -L$(ACCFFT_LIB) $(ACCFFT_LIBNAME) -Wl,-rpath,$(ACCFFT_LIB)

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
# OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
# OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)

################################################################################

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_DPREC_A2A) $(TARGET_EXE_DPREC_NB) $(TARGET_EXE_ISR)

all: default

$(TARGET_EXE_ISR): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_ISR)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_A2A)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_NB)
	
$(TARGET_EXE_DPREC_A2A): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_DPREC_A2A)

$(TARGET_EXE_DPREC_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB) $(LIB_DPREC_NB)

install: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_ISR)
	mkdir -p $(PREFIX)/bin
	cp $(TARGET_EXE_ISR) $(PREFIX)/bin
	cp $(TARGET_EXE_A2A) $(PREFIX)/bin
	cp $(TARGET_EXE_NB) $(PREFIX)/bin

test:
	@echo $(SRC)
	@echo $(INC)
	@echo $(LIB)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_ISR)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	
destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET_EXE_ISR)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)
	rm -f $(OBJ_DIR)/*
	
info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
