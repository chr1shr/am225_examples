# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=diffuse.o ho_transport.o
src=$(patsubst %.o,%.cc,$(objs))
execs=d_solve d_cons_test ho_solve

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
	$(cxx) $(cflags) -c $< $(ode_iflags)

ho_solve: ho_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

d_solve: d_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

d_cons_test: d_cons_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend
