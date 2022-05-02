# Prevent the parallel execution when calling this Makefile
.NOTPARALLEL:

.PHONY: info
info: logo
	$(info prefix = $(PREFIX)/lib )
	$(info compiler = $(shell $(CXX) --version))
	$(info compil. flags = $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD)
	$(info linker flags = -shared $(LDFLAGS))
	$(info compil. options = $(OPTS))
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
	$(info H3LPR:)
	$(info - include: -I$(H3LPR_INC) )
	$(info - lib: -L$(H3LPR_LIB) $(H3LPR_LIBNAME) -Wl,-rpath,$(H3LPR_LIB))
	$(info ------------)
	$(info LIST OF OBJECTS:)
	$(info - SRC = $(SRC))
	$(info - OBJ A2A = $(OBJ_A2A))
	$(info - OBJ NB = $(OBJ_NB))
	$(info - DEP = $(DEP))
	$(info - LGF_DATA = $(LGF_DATA))
	$(info ------------)

.NOTPARALLEL: logo
.PHONY: logo
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
	@echo "    	(C) UCLouvain - Appache 2.0                    "
	@echo "----------------------------------------------------"