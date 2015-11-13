###############################################
## Make file for Heston model discretization ##
###############################################


##### Insert here location of headers and libraries

# headers
ARMADILLO_INC	= /media/mircea/Stuff/Mircea/Programming/C++/libraries/armadillo-6.100.1/include

# libraries
OPENBLAS_LIB 	= /usr/lib/openblas/lib
LAPACK_LIB		= /usr/lib/lapack

#################################################################

CC				= g++
DEBUG			= -DARMA_NO_DEBUG -DNDEBUG
PROG            = heston
OBJS            = heston.o
FILES			= bpm.cpp
INCLUDES        = -I $(ARMADILLO_INC)
LIBS            = -L $(OPENBLAS_LIB) -lopenblas -L $(LAPACK_LIB) -llapack
CFLAGS        	= -O3 -std=c++11 -fopenmp -march=native $(DEBUG)

build: $(FILES)
	$(CC) -o $(PROG) $(FILES) $(INCLUDES) $(LIBS) $(CFLAGS)

clean:
	rm -f *.o core
