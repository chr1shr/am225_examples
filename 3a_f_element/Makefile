# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=conj_grad.o cubic_1d_fe.o quadrat.o
src=$(patsubst %.o,%.cc,$(objs))
execs=fe_1d_test fe_1d_conv q2d_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

fe_1d_test: fe_1d_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

fe_1d_conv: fe_1d_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

q2d_test: q2d_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend
