CC  = g++

CURPATH := $(shell pwd)
SUBDIRS := $(shell find $(CURPATH) -maxdepth 1 -type d)
PRJDIRS := $(basename $(patsubst $(CURPATH)/%, %, $(SUBDIRS)))
# PRJDIRS := $(patsubst $(CURPATH)/%, %, $(SUBDIRS))
PRJDIRS := $(filter-out $(CURPATH), $(PRJDIRS))
PRJDIRS := $(filter-out .vscode, $(PRJDIRS))
PRJDIRS := $(filter-out test, $(PRJDIRS))
PRJDIRS := $(filter-out cpp-api, $(PRJDIRS))
PRJDIRS := $(filter-out matlab-dll, $(PRJDIRS))
PRJDIRS := $(filter-out matlab-mex, $(PRJDIRS))

MKLROOT ?= /opt/intel/oneapi/mkl/latest
EIGROOT ?= /opt/eigen-3.4.0
BLAROOT = ./skewblas

LIBDIRS = $(MKLROOT)/lib

INCDIRS = -I $(CURPATH) -I $(MKLROOT)/include -I $(BLAROOT) -I $(EIGROOT)

LIBS = -Wl,--start-group $(LIBDIRS)/libmkl_intel_ilp64.a \
                  $(LIBDIRS)/libmkl_sequential.a \
                  $(LIBDIRS)/libmkl_core.a \
       -Wl,--end-group

EXTRA = -O3
# EXTRA = -g

LFLAGS = $(EXTRA) -m64 $(LIBS) -lpthread -lm -ldl
CFLAGS = $(EXTRA) -DMKL_ILP64 -m64 $(INCDIRS) -c

target = dexp

# ----------------------------------------------------------------------
# Core sources in the top directory
# ----------------------------------------------------------------------
SRC      = $(wildcard *.cpp)
OBJ      = $(SRC:.cpp=.o)

# Exclude main.cpp when linking test executables (if main.cpp exists)
CORE_SRC = $(filter-out main.cpp,$(SRC))
CORE_OBJ = $(CORE_SRC:.cpp=.o)

# ----------------------------------------------------------------------
# Subproject sources (unchanged)
# ----------------------------------------------------------------------
SUBSRC := $(wildcard $(addsuffix /*.cpp, $(PRJDIRS)))
SUBOBJ := $(patsubst %.cpp,%.o,$(SUBSRC))

export

# ----------------------------------------------------------------------
# Test infrastructure: one executable per test/*.cpp
# ----------------------------------------------------------------------
TESTDIR  := test
TESTSRC  := $(wildcard $(TESTDIR)/*.cpp)
# Test names = basename without directory and extension
TESTS    := $(patsubst $(TESTDIR)/%.cpp,%,$(TESTSRC))

.PHONY: all lib rebuild remake clean tests $(TESTS)

all : $(target)

# Generic object rule (top-level sources)
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

# Dependency files for top-level sources
%.d: %.cpp
	@set -e; rm -f $@; $(CC) $(INCDIRS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

sinclude $(SRC:.cpp=.d)

# ----------------------------------------------------------------------
# Main executable: uses main.cpp + all other objects (as before)
# ----------------------------------------------------------------------
$(target) : $(OBJ)
	@for dir in ${PRJDIRS}; do $(MAKE) -C $(CURPATH)/$$dir || exit; done
	$(CC) $^ $(SUBOBJ) -o $@ $(LFLAGS)

# ----------------------------------------------------------------------
# Build all libraries in subprojects (no executable)
# ----------------------------------------------------------------------
lib :
	@for dir in ${PRJDIRS}; do $(MAKE) -C $(CURPATH)/$$dir || exit; done

# ----------------------------------------------------------------------
# Rebuild the main executable
# ----------------------------------------------------------------------
rebuild : $(SRC)
	rm -rf $(target)
	@for dir in ${PRJDIRS}; do $(MAKE) -C $(CURPATH)/$$dir || exit; done
	$(CC) $(OBJ) $(SUBOBJ) -o $(target) $(LFLAGS)

remake: clean
	@echo "Rebuilding project..."
	$(MAKE)

# ----------------------------------------------------------------------
# Tests: each test/*.cpp -> executable named after the source file
#
#   test/foo.cpp  ->  executable ./foo
#
# They link against CORE_OBJ (no main.cpp) + SUBOBJ, so main.cpp is
# NEVER overridden or duplicated.
# ----------------------------------------------------------------------
tests: $(TESTS)

# Object files for tests live under test/
$(TESTDIR)/%.o: $(TESTDIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

# Each test target: name "foo" builds from test/foo.o + CORE_OBJ + SUBOBJ
$(TESTS): %: $(TESTDIR)/%.o $(CORE_OBJ)
	@for dir in ${PRJDIRS}; do $(MAKE) -C $(CURPATH)/$$dir || exit; done
	$(CC) $^ $(SUBOBJ) -o $@ $(LFLAGS)

# ----------------------------------------------------------------------
# Clean
# ----------------------------------------------------------------------
clean :
	@for dir in ${PRJDIRS}; do $(MAKE) -C $(CURPATH)/$$dir clean || exit; done
	rm -rf $(target) *.o *.d *.d.* $(TESTDIR)/*.o $(TESTS)
