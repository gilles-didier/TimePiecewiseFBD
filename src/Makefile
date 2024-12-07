# list of executable files to produce
GET_ML = getML
# .o files necessary to build the executables


OBJ_GET_ML =  MyR.o MyRandom.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MinimizeNLOpt.o MCMCImportanceSamplingFossilInt.o Distribution.o FBDDensity.o GetML.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = GET_ML

# compiler and its options
CC = gcc
#CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT -g 
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -I/usr/include/gsl -I/usr/share/R/include -g -fopenmp

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(GET_ML)

# build the executable
$(GET_ML): $(OBJ_GET_ML)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -lnlopt -fopenmp  -g


############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release
