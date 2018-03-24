# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=pss_read.o
src=$(patsubst %.o,%.cc,$(objs))
execs=pss_test

all: $(objs) $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $<

pss_test: pss_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^

.PHONY: clean depend
