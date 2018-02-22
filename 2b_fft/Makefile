# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=poisson_fft.o file_output.o sol_base.o sol_rk4.o kdv.o
src=$(patsubst %.o,%.cc,$(objs))
execs=ranal pfft_test pfft_conv fft_deriv kdv_test

all: $(objs) $(execs)

# Include the file dependencies
#-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $<

ranal: ranal.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

pfft_test: pfft_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

pfft_conv: pfft_conv.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

fft_deriv: fft_deriv.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

kdv_test: kdv_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags)

.PHONY: clean depend
