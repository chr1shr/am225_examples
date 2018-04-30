# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=dop853.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=d_test brus_conv

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

d_test: d_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

brus_conv: brus_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
