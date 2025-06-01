#CC  = g++

#MKLROOT = /opt/intel/oneapi/mkl/latest

#LIBDIRS = $(MKLROOT)/lib

#INCDIRS = -I $(MKLROOT)/include

#LIBS = -Wl,--start-group $(LIBDIRS)/libmkl_intel_ilp64.a $(LIBDIRS)/libmkl_sequential.a $(LIBDIRS)/libmkl_core.a -Wl,--end-group

#LFLAGS = $(EXTRA) -m64 $(LIBS) -lpthread -lm -ldl
#CFLAGS = $(EXTRA) -DMKL_ILP64  -m64 $(INCDIRS) -c

# EXTRA = -g
#EXTRA = -O3



# target = skewsymm

SRC = $(wildcard *.cpp)

OBJ = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all : $(OBJ) 


%.o:%.cpp
	$(CC) $(CFLAGS) $< -o $@   

%.d:%.cpp
	@set -e; rm -f $@; $(CC) $(INCDIRS) -MM $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

sinclude $(SRC:.cpp=.d)

# $(target) : $(OBJ)
# 	$(CC) $^ -o $@  $(LFLAGS)

.PHONY : clean
clean :
	rm -rf $(target) *.o *.d *.d.*
