# ================================================================
# Makefile for libskewblas.a
# ================================================================

CC       = g++
EXTRA    = -O3
CFLAGS   = $(EXTRA) -DMKL_ILP64 -m64 -fPIC
AR       = ar rcs

MKLROOT  = /opt/intel/oneapi/mkl/latest
EIGROOT  = /opt/eigen-3.4.0

INCDIRS  = -I$(MKLROOT)/include -I$(EIGROOT) -I.
SRC      = $(wildcard *.cpp)
OBJ      = $(SRC:.cpp=.o)
TARGET   = libskewblas.a

all: $(TARGET)

$(TARGET): $(OBJ)
	$(AR) $@ $(OBJ)
	@echo "Created static library: $(TARGET)"

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCDIRS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
	@echo "Cleaned libskewblas.a"

.PHONY: all clean



# # target = skewsymm

# SRC = $(wildcard *.cpp)

# OBJ = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

# all : $(OBJ) 


# %.o:%.cpp
# 	$(CC) $(CFLAGS) $< -o $@   

# %.d:%.cpp
# 	@set -e; rm -f $@; $(CC) $(INCDIRS) -MM $< > $@.$$$$; \
# 	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
# 	rm -f $@.$$$$

# sinclude $(SRC:.cpp=.d)

# # $(target) : $(OBJ)
# # 	$(CC) $^ -o $@  $(LFLAGS)

# .PHONY : clean
# clean :
# 	rm -rf $(target) *.o *.d *.d.*
