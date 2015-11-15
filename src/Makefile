###############################################
## Make file for Heston model discretization ##
###############################################


##### Insert here location of headers and libraries
###################################################

# headers
ARMADILLO_INC	= /media/mircea/Stuff/Mircea/Programming/C++/libraries/armadillo-6.100.1/include

# libraries
OPENBLAS_LIB 	= /usr/lib/openblas/lib
LAPACK_LIB		= /usr/lib/lapack

###################################################
###################################################

CC				= g++
DEBUG			= -DARMA_NO_DEBUG -DNDEBUG
PROG            = heston-mc
OBJECTS			= heston-mc.o hestonDiscr.o hestonEuler.o hestonMilstein.o KahlJackel.o option.o payoff.o transformedVolatility.o
FILES			= 
INCLUDES        = -I $(ARMADILLO_INC)
LIBS            = -L $(OPENBLAS_LIB) -lopenblas -L $(LAPACK_LIB) -llapack
CFLAGS        	= -O3 -std=c++11 -fopenmp -march=native $(DEBUG)

heston-mc: $(OBJECTS)
	$(CC) -o $(PROG) $(OBJECTS) $(INCLUDES) $(LIBS) $(CFLAGS)


## Generic rules:
.SUFFIXES:          # get rid of that annoying Modula rule
.SUFFIXES: .cpp .h .o

#
# Rules
#
.cpp.o:
	$(CC) -c $*.cpp $(INCLUDES) $(LIBS) $(CFLAGS)
#
# Dependencies
#
heston-mc.o: option.h payoff.h helpers.h hestonDiscr.h hestonEuler.h hestonMilstein.h KahlJackel.h transformedVolatility.h
#
payoff.o: payoff.h
#
option.o: option.h payoff.h
#
hestonDiscr.o: hestonDiscr.h option.h
#
hestonEuler.o: hestonEuler.h hestonDiscr.h option.h
#
hestonMilstein.o: hestonMilstein.h hestonDiscr.h option.h
#
KahlJackel.o: KahlJackel.h hestonDiscr.h option.h
#
transformedVolatility.o: transformedVolatility.h hestonDiscr.h option.h


clean:
	\rm *.o *~ heston-mc
