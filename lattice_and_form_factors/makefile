# MAKEFILE for lattice and form factors unit "tests"

#--------------------------------------General------------------------------------------
CXX 	:= g++ # This is the main compiler
#CXX 	:= icc # This is the main compiler
SRCDIR 	:= src
HEADDIR := include
BUILDDIR := build
TARGET 	:= bin/run


#--------------------------------------Sources and header files------------------------------------------
SRCEXT 	:= cpp
HEADEXT := h
#SOURCES := $(shell find $(SRCDIR) -name '*.$(SRCEXT)') 
#HEADERS := $(shell find $(HEADDIR) -name '*.$(HEADEXT)') 
SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT)) 
HEADERS := $(shell ls -1 $(HEADDIR)/*.$(HEADEXT)) 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))


#--------------------------------------Compiler settings------------------------------------------

CFLAGS += -std=c++17  -g  # 				General compiler flags
#CFLAGS += -D BOOST_DISABLE_ASSERTS # 				Disable boost assert checks
DBFLAGS := -O2 -g # 						Compiler flags for debugging
PROFFLAGS := -O2 -g # 						Compiler flags for profiling
RUNFLAGS := -O3 # 				Compiler flags for quick compile and run
OPTIMIZEFLAGS := -flto -march=native -O3  # Compiler flags for optimal speed
#OPTIMIZEFLAGS := -xHOST -fp-model fast -lpthread -O3 -openmp -DOPENMP # -opt-prefetch 0 # Compiler flags for optimal speed
LIB := -lgsl -lgslcblas -lboost_system -lboost_timer -lm #				Library flags
INC := -I include #						Additional include paths


#--------------------------------------Targets ------------------------------------------
$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@mkdir -p bin
	@mkdir -p dat
	@echo " $(CXX) $^ -o $(TARGET) $(LIB)"; $(CXX) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

run: 	CFLAGS += $(RUNFLAGS)
run: 	$(TARGET)

debug: 	CFLAGS += $(DBFLAGS)
debug: 	$(TARGET)

prof: 	CFLAGS += $(PROFFLAGS)
prof: 	$(TARGET)

optimize: CFLAGS += $(OPTIMIZEFLAGS)
optimize: $(SOURCES) $(HEADERS)
	@mkdir -p bin
	@mkdir -p dat
	$(CXX) $(CFLAGS) $(INC) $(SOURCES) -o $(TARGET) $(LIB)
	#rm *.o

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean

