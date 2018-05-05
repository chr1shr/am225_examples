# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=sol_base.o sol_euler.o sol_improv_e.o sol_imp_onestep.o \
	 sol_sep.o orbit.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=orb_solve orb_conv

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

orb_solve: orb_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

orb_conv: orb_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
