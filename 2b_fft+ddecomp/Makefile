# Include the configuration file
include ../config.mk

cflags+=-std=c++11

# Lists of files to be built
objs=poisson_fft.o schwarz.o conj_grad.o schur.o file_output.o sol_base.o sol_rk4.o kdv.o
src=$(patsubst %.o,%.cc,$(objs))
execs=music_fft fftw_time pfft_test pfft_conv fft_deriv kdv_test sch_test schur_test

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

music_fft: music_fft.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

pfft_test: pfft_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

pfft_conv: pfft_conv.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

fftw_time: fftw_time.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

fft_deriv: fft_deriv.cc
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

kdv_test: kdv_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

sch_test: sch_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

schur_test: schur_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

.PHONY: clean depend
