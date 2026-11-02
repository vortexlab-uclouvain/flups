################################################################################
# @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
#      See LICENSE file in top-level directory
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_clang
include ../../$(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
NAME := flups
# executable naming
TARGET_EXE_ISR := $(NAME)_validation_isr
TARGET_EXE_A2A := $(NAME)_validation_a2a
TARGET_EXE_NB := $(NAME)_validation_nb

#-----------------------------------------------------------------------------
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FLUPS
FLUPS_INC ?= ../../include
FLUPS_LIB ?= ../../lib
INC += -I$(FLUPS_INC)

#---- MPI
MPI_DIR ?= /usr
MPI_INC ?= $(MPI_DIR)/include
MPI_LIB ?= $(MPI_DIR)/lib
INC += -I$(MPI_INC)
LIB += -L$(MPI_LIB) -Wl,-rpath,$(MPI_LIB)

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
	ifneq (,$(findstring nvcc,$(shell $(CXX) --version)))
		CXXFLAGS += --forward-unknown-to-host-compiler
		LDFLAGS  += --forward-unknown-to-host-compiler
	endif
endif

#---- Use HIP
ifneq (,$(findstring -DBACKEND_HIP,$(OPTS)))
	HIP_INC ?= /opt/rocm/include
	HIP_LIB ?= /opt/rocm/lib
	HIP_LIBNAME ?= -lrocblas -lhipfft -lamdhip64
	INC += -I$(HIP_INC)
	LIB += -L$(HIP_LIB) $(HIP_LIBNAME) -Wl,-rpath,$(HIP_LIB)
	CXXFLAGS += -x hip
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

$(TARGET_EXE_ISR): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_isr -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_A2A): $(OBJ)
	$(CXX) $(LDFLAGS)  $^ -o $@ -L$(FLUPS_LIB) -lflups_a2a -Wl,-rpath,$(FLUPS_LIB) $(LIB)

$(TARGET_EXE_NB): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ -L$(FLUPS_LIB) -lflups_nb -Wl,-rpath,$(FLUPS_LIB) $(LIB)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE_ISR)
	rm -f $(TARGET_EXE_A2A)
	rm -f $(TARGET_EXE_NB)

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
