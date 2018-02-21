# Include the configuration file
cxx=g++-mp-7
cflags=-fopenmp -Wall -ansi -pedantic

# Lists of files to be built
objs=dop853.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=d_test brus_conv

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

d_test: d_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

brus_conv: brus_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
