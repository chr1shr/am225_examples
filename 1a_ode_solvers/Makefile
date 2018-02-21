# Include the configuration file
cxx=g++-mp-7
cflags=-fopenmp -Wall -ansi -pedantic

# Lists of files to be built
objs=sol_base.o sol_euler.o sol_ralston.o sol_heun3.o sol_rk4.o sol_hammer_h.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=ode_test ode_conv

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

ode_test: ode_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

ode_conv: ode_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
