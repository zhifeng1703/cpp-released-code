CC  = g++

CURPATH := $(shell pwd)
SUBDIRS := $(shell find $(CURPATH) -maxdepth 1 -type d)
PRJDIRS := $(basename $(patsubst $(CURPATH)/%, %, $(SUBDIRS)))
# PRJDIRS := $(patsubst $(CURPATH)/%, %, $(SUBDIRS))
PRJDIRS := $(filter-out $(CURPATH), $(PRJDIRS))
PRJDIRS := $(filter-out .vscode, $(PRJDIRS))
PRJDIRS := $(filter-out test, $(PRJDIRS))


MKLROOT = /opt/intel/oneapi/mkl/latest
BLAROOT = ./blas_lapack_interface

LIBDIRS = $(MKLROOT)/lib

INCDIRS = -I $(MKLROOT)/include -I $(BLAROOT)

LIBS = -Wl,--start-group $(LIBDIRS)/libmkl_intel_ilp64.a $(LIBDIRS)/libmkl_sequential.a $(LIBDIRS)/libmkl_core.a -Wl,--end-group


EXTRA = -O3
#EXTRA = -g


LFLAGS = $(EXTRA) -m64 $(LIBS) -lpthread -lm -ldl
CFLAGS = $(EXTRA) -DMKL_ILP64  -m64 $(INCDIRS) -c


target = skewsymm



SRC = $(wildcard *.cpp)

OBJ := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

SUBSRC := $(wildcard $(addsuffix *.cpp, $(PRJDIRS)/))
SUBOBJ := $(patsubst %.cpp,%.o,$(SUBSRC))

export

all : $(target)

%.o:%.cpp
	$(CC) $(CFLAGS) $< -o $@   

%.d:%.cpp
	@set -e; rm -f $@; $(CC) $(INCDIRS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

sinclude $(SRC:.cpp=.d)

$(target) : $(OBJ)
	@for dir in ${PRJDIRS}; do make -C $(CURPATH)/$$dir ||exit; done
	$(CC) $^ $(SUBOBJ) -o $@  $(LFLAGS)

lib :
	@for dir in ${PRJDIRS}; do make -C $(CURPATH)/$$dir ||exit; done

rebuild : $(SRC)
	rm -rf $(target)
	@for dir in ${PRJDIRS}; do make -C $(CURPATH)/$$dir ||exit; done
	$(CC) $(OBJ) $(SUBOBJ) -o $(target)  $(LFLAGS)

.PHONY : remake
remake : 
	make clean
	make


.PHONY : clean
clean :
	@for dir in ${PRJDIRS}; do make -C $(CURPATH)/$$dir clean ||exit; done
	rm -rf $(target) *.o *.d *.d.*
