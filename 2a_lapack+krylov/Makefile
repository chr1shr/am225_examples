# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=file_output.o conj_grad.o poisson.o rbf.o lp_solve.o
src=$(patsubst %.o,%.cc,$(objs))
execs=blas_test blas_test2 rbf_test rbf_test2 basic_cg_test rbf_time

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

liblinal.a: $(objs)
	rm -f $@
	ar rs $@ $^

%.o: %.cc
	$(cxx) $(cflags) -c $<

blas_test: blas_test.cc
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

blas_test2: blas_test2.cc
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

basic_cg_test: basic_cg_test.cc liblinal.a
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

rbf_test: rbf_test.cc liblinal.a
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

rbf_test2: rbf_test2.cc liblinal.a
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

rbf_time: rbf_time.cc liblinal.a
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend
