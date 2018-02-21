# Include the configuration file
cxx=g++-mp-7
cflags=-fopenmp -Wall -ansi -pedantic

# Lists of files to be built
objs=sol_rk4d.o sol_extrap.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=dens_test extrap_test extrap_conv extrap_ctri

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

dens_test: dens_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

extrap_conv: extrap_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

extrap_ctri: extrap_ctri.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

extrap_test: extrap_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
