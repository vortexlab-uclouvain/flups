################################################################################
# @copyright Copyright (c) Université catholique de Louvain (UCLouvain), Belgique 
#      See LICENSE file in top-level directory
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_clang

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-------------------------------------------------------------------------------
# make predefines CXX (g++), so a plain ?= would never apply
ifeq ($(origin CXX),default)
CXX := mpic++
endif
AR ?= ar

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
# git commit
GIT_COMMIT ?= $(shell git rev-parse --short HEAD)

PREFIX ?= ./
NAME := flups

#-------------------------------------------------------------------------------
SRC_DIR := ./src
OBJ_DIR := ./build
KERNEL_DIR := ./kernel

# library naming
TARGET_LIB_ISR := $(OBJ_DIR)/lib$(NAME)_isr
TARGET_LIB_A2A := $(OBJ_DIR)/lib$(NAME)_a2a
TARGET_LIB_NB  := $(OBJ_DIR)/lib$(NAME)_nb

## add the headers to the includes
INC := -I$(SRC_DIR)

#-------------------------------------------------------------------------------
# Required Dependencies
#---- HDF5
HDF5_DIR ?= /usr
HDF5_INC ?= $(HDF5_DIR)/include
HDF5_LIB ?= $(HDF5_DIR)/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) $(RPATH_FLAG)$(HDF5_LIB)$(RPATH_END)

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
# Backend -- FFTW by default
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
CXXFLAGS += -x cu
ifneq (,$(IS_NVCC))
CXXFLAGS += --forward-unknown-to-host-compiler --extended-lambda --expt-relaxed-constexpr
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
CXXFLAGS += -x hip
endif

#-------------------------------------------------------------------------------
# LGF kernels: by default the kernel data is installed in the include directory
LGF_PATH := $(abspath $(PREFIX)/include)
DEF += -DKERNEL_PATH=\"$(LGF_PATH)\"
LGF_DATA := $(wildcard $(KERNEL_DIR)/*.ker)

#-------------------------------------------------------------------------------
## sources and headers
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)
API := $(wildcard $(SRC_DIR)/*.h)

## generate object list, one set per communication scheme
OBJ_ISR := $(SRC:%.cpp=$(OBJ_DIR)/isr_%.o)
OBJ_A2A := $(SRC:%.cpp=$(OBJ_DIR)/a2a_%.o)
OBJ_NB := $(SRC:%.cpp=$(OBJ_DIR)/nb_%.o)
DEP := $(OBJ_ISR:%.o=%.d) $(OBJ_A2A:%.o=%.d) $(OBJ_NB:%.o=%.d)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)

################################################################################
# mandatory flags
M_FLAGS := -fPIC -DGIT_COMMIT=\"$(GIT_COMMIT)\"

#-------------------------------------------------------------------------------
$(OBJ_DIR)/isr_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_ISR $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/nb_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) -DCOMM_NONBLOCK $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/a2a_%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -E $< -o $@

################################################################################
.PHONY: default
default:
	@$(MAKE) info
	@$(MAKE) lib_static

# for the validation and the tests, do a static lib
.PHONY: validation test python
validation: install_static
test: install_static
python: install_dynamic

#-------------------------------------------------------------------------------
# libraries
.PHONY: all all2all nonblocking lib lib_static lib_dynamic preproc
all: lib_static lib_dynamic
all2all: $(TARGET_LIB_A2A).a $(TARGET_LIB_A2A).so
nonblocking: $(TARGET_LIB_NB).a $(TARGET_LIB_NB).so $(TARGET_LIB_ISR).a $(TARGET_LIB_ISR).so
lib: lib_static
lib_static: $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a $(TARGET_LIB_ISR).a
lib_dynamic: $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so $(TARGET_LIB_ISR).so
preproc: $(IN)

$(TARGET_LIB_ISR).so: $(OBJ_ISR)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_A2A).so: $(OBJ_A2A)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_NB).so: $(OBJ_NB)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB_ISR).a: $(OBJ_ISR)
	$(AR) rvs $@ $^

$(TARGET_LIB_A2A).a: $(OBJ_A2A)
	$(AR) rvs $@ $^

$(TARGET_LIB_NB).a: $(OBJ_NB)
	$(AR) rvs $@ $^

#-------------------------------------------------------------------------------
# installation
.PHONY: install install_static install_dynamic install_headers install_macos lib_macos
install:
	@$(MAKE) info
	@$(MAKE) install_static

install_headers:
	@mkdir -p $(PREFIX)/lib $(PREFIX)/include
	@cp $(API) $(LGF_DATA) $(PREFIX)/include

install_static: lib_static install_headers
	@cp $(TARGET_LIB_ISR).a $(TARGET_LIB_A2A).a $(TARGET_LIB_NB).a $(PREFIX)/lib

install_dynamic: lib_dynamic install_headers
	@cp $(TARGET_LIB_ISR).so $(TARGET_LIB_A2A).so $(TARGET_LIB_NB).so $(PREFIX)/lib

# macOS specific: create a unified dynamic library from the static all-to-all library
lib_macos: $(TARGET_LIB_A2A).a
	@mkdir -p $(PREFIX)/lib
	@echo "Creating unified macOS dynamic library..."
	$(CXX) -dynamiclib -o $(PREFIX)/lib/libflups.dylib \
		-Wl,-force_load,$(TARGET_LIB_A2A).a \
		$(LIB) $(LDFLAGS) \
		-install_name @rpath/libflups.dylib

install_macos: install_static lib_macos
	@echo "macOS unified dynamic library created: $(PREFIX)/lib/libflups.dylib"

#-------------------------------------------------------------------------------
.PHONY: clean destroy
clean:
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TARGET_LIB_ISR).so $(TARGET_LIB_ISR).a
	@rm -f $(TARGET_LIB_A2A).so $(TARGET_LIB_A2A).a
	@rm -f $(TARGET_LIB_NB).so $(TARGET_LIB_NB).a

destroy: clean
	@rm -f $(OBJ_DIR)/*.d $(OBJ_DIR)/*.in

#-------------------------------------------------------------------------------
# mentioning this target will export all the current variables to child-make processes
.EXPORT_ALL_VARIABLES:

.PHONY: info
info:
	@$(MAKE) --file=info.mak

#-------------------------------------------------------------------------------
-include $(DEP)

# mkdir the needed dirs
$(shell mkdir -p $(OBJ_DIR))
