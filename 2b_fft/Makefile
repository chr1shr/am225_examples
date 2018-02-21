# Include the configuration file
cxx=g++-mp-7
cflags=-fopenmp -O3 -Wall -ansi -pedantic
lflags=-L/opt/local/lib -lfftw3

# Lists of files to be built
#objs=file_output.o conj_grad.o poisson.o rbf.o lp_solve.o
#src=$(patsubst %.o,%.cc,$(objs))
execs=ranal

all: $(execs)

# Include the file dependencies
#-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
#clean:
#	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

ranal: ranal.cc
	$(cxx) $(cflags) -o $@ $^ $(lflags)

.PHONY: clean depend
