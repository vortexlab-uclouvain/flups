################################################################################
# ARCH DEPENDENT VARIABLES
include make_arch/make.vagrant_intel

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
NAME := flups
TARGET_EXE := $(NAME)_validation
TARGET_LIB := build/lib$(NAME).so

PREFIX ?= ./

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

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD -c $< -o $@

################################################################################
default: $(TARGET_EXE)

all: $(TARGET_EXE) $(TARGET_LIB)

lib: $(TARGET_LIB)

$(TARGET_EXE): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

$(TARGET_LIB): $(OBJ)
	$(CXX) -shared $(LDFLAGS) $^ -o $@ $(LIB)

install: $(TARGET_LIB)
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	cp $(TARGET_LIB) $(PREFIX)/lib
	cp $(HEAD) $(PREFIX)/include

test:
	@echo $(SRC)

clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_LIB)

destroy:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.d
	rm -f $(TARGET_EXE)
	rm -f $(TARGET_LIB)
	rm -f $(OBJ_DIR)/*
	rm -f include/*
	rm -f lib/*

info:
	$(info SRC = $(SRC))
	$(info OBJ = $(OBJ))
	$(info DEP = $(DEP))

-include $(DEP)
