# Include the configuration file
include ../config.mk

# Lists of files to be built
#objs=file_output.o conj_grad.o poisson.o rbf.o lp_solve.o
#src=$(patsubst %.o,%.cc,$(objs))
execs=ranal

all: $(execs)

# Include the file dependencies
#-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $<

ranal: ranal.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

.PHONY: clean depend
