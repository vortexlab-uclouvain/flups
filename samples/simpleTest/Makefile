################################################################################
# @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
#      See LICENSE file in top-level directory
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_gcc

include ../../$(ARCH_FILE)

FLUPSDIR := ../../

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := simple_test
# executable naming
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
FLUPS_INC ?= ../../include
FLUPS_LIB ?= ../../lib
INC += -I$(FLUPS_INC)
# LIB += -L$(FLUPS_LIB) $(FLUPS_LIBNAME) -Wl,-rpath,$(FLUPS_LIB)

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


#---- H3LPR
H3LPR_INC ?= /usr/include
H3LPR_LIB ?= /usr/lib
H3LPR_LIBNAME ?= -lh3lpr
INC += -I$(H3LPR_INC)
LIB += -L$(H3LPR_LIB) $(H3LPR_LIBNAME) -Wl,-rpath,$(H3LPR_LIB)


#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
# HEAD := $(wildcard $(SRC_DIR)/*.hpp)

## generate object list
# DEP := $(SRC:%.c=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

################################################################################

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE_A2A) $(TARGET_EXE_NB)

all: $(TARGET_EXE_A2A) $(TARGET_EXE_NB) $(TARGET_EXE_DPREC_A2A) $(TARGET_EXE_DPREC_NB)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_DPREC_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_dprec_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_DPREC_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_dprec_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB)

install: $(TARGET_EXE_A2A) 
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
	rm -f $(TARGET_EXE_DPREC_A2A)
	rm -f $(TARGET_EXE_DPREC_NB)
	
destroy: clean 
	rm -f $(OBJ_DIR)/*
	
info:
	@echo $(ARCH_FILE)
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	# $(info DEP = $(DEP))

# -include $(DEP)
