# MAKEFILE for BosonFlow

#--------------------------------------General------------------------------------------
CXX 	:= g++
SRCDIR  := src
SRCDIRLATTICE := lattice_and_form_factors/src
HEADDIRLATTICE := lattice_and_form_factors/include
SRCDIRCUBATURE := cubature
HEADDIRCUBATURE := cubature
HEADDIR := include
BUILDDIR := build
TARGET := ./bin/run

# Files to exclude during compilation
EXCLUDE_CUBATURE := cubature/test.c cubature/clencurt_gen.c

#--------------------------------------Sources and header files------------------------------------------
SRCEXT  := cpp
CSRCEXT := c
HEADEXT := h
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
SOURCES_LATTICE := $(wildcard $(SRCDIRLATTICE)/*.$(SRCEXT))
SOURCES_CUBATURE := $(wildcard $(SRCDIRCUBATURE)/*.$(CSRCEXT))
SOURCES_CUBATURE := $(filter-out $(EXCLUDE_CUBATURE), $(SOURCES_CUBATURE))
HEADERS := $(shell find $(HEADDIR) $(HEADDIRLATTICE) $(HEADDIRCUBATURE) -type f -name '*.$(HEADEXT)')
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS2 := $(patsubst $(SRCDIRLATTICE)/%,$(BUILDDIR)/%,$(SOURCES_LATTICE:.$(SRCEXT)=.o))
OBJECTS3 := $(patsubst $(SRCDIRCUBATURE)/%,$(BUILDDIR)/%,$(SOURCES_CUBATURE:.$(CSRCEXT)=.o))
ALL_OBJECTS := $(OBJECTS) $(OBJECTS2) $(OBJECTS3)

include config.mk

#CFLAGS += -D USE_KNN_SEARCH

# -- Julia integration --------------------------------------------------------
ifdef JULIA_DIR_LA
CFLAGS += -D JULIA_AVAILABLE
endif

#--------------------------------------Compiler settings------------------------------------------

# compiler settings for cluster in tue
CFLAGS += -std=c++20 -D OPENMP -fopenmp #			General compiler flags	
# Eigen stack allocation limit in bytes. 512 bytes corresponds to up to about 100 float temporaries per thread.
CFLAGS += -DEIGEN_STACK_ALLOCATION_LIMIT=512
DBFLAGS := -g #							Compiler flags for debugging
RUNFLAGS := -D BOOST_DISABLE_ASSERTS  -DNDEBUG \
            -D BOOST_BIND_GLOBAL_PLACEHOLDERS \
            -D BOOST_ALLOW_DEPRECATED_HEADERS #			Compiler flags for quick compile and run
OPTIMIZEFLAGSCLANG := -O3 #-Ofast (alternative)
OPTIMIZEFLAGSGPP :=   -march=native -O3 #-flto -march=native -O3
LIB := -march=native -Wdeprecated-declarations -Wl,-rpath=$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib #-flto -march=native -Wdeprecated-declarations -Wl,-rpath=$(HDF5_DIR)/lib -L$(HDF5_DIR)/lib

ifdef JULIA_DIR_LA 
	LIB += -L$(JULIA_DIR)/lib -Wl,-rpath,$(JULIA_DIR)/lib 
endif
LIB += -lhdf5 -lm -fopenmp -lfftw3 -fPIC
ifdef JULIA_DIR_LA
	LIB += -ljulia # Library flags  # for IR basis integration
endif						
INC += -I $(HEADDIR)  -I $(HEADDIRLATTICE)  -I $(HEADDIRCUBATURE)
ifdef JULIA_DIR_LA
INC += -I$(JULIA_DIR)/include/julia			#Additional include paths
endif
#--------------------------------------Targets ------------------------------------------
$(TARGET): $(ALL_OBJECTS)
	@echo " Linking..."
	@mkdir -p bin
	@mkdir -p dat
	@echo " $(CXX) $^ -o $(TARGET) $(LIB)"; $(CXX) $^ -o $(TARGET) $(LIB)

$(OBJECTS2): $(BUILDDIR)/%.o: $(SRCDIRLATTICE)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(dir $@)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

$(OBJECTS3): $(BUILDDIR)/%.o: $(SRCDIRCUBATURE)/%.$(CSRCEXT) $(HEADERS)
	@mkdir -p $(dir $@)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

$(OBJECTS): $(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(dir $@)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

run: 	CFLAGS += $(OPTIMIZEFLAGSGPP) $(RUNFLAGS) -msse
run: 	$(TARGET)

runclang: CXX = clang++ 
runclang: CFLAGS += $(OPTIMIZEFLAGSCLANG) $(RUNFLAGS) -msse -w
runclang: $(TARGET)

dbg:    CFLAGS += -D DEBUG_BUILD $(DBFLAGS) -msse
dbg:    $(TARGET)

dbgclang:    CXX = clang++
dbgclang:    CFLAGS += -D DEBUG_BUILD $(DBFLAGS) -msse
dbgclang:    $(TARGET)


clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: run runclang dbg dbgclang clean

