# Load the common configuration file
include ../../config.mk

iflags=-I../../tgmg -I..

objs=mgs_fem_vr.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mgs_test

all:
	$(MAKE) -C ../../tgmg
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

mgs_test: mgs_test.cc mgs_fem_vr.o
	$(cxx) $(cflags) $(iflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs)

clean-all:
	$(MAKE) -C ../../tgmg clean
	rm -f $(execs) $(objs)

.PHONY: clean all clean-all executables depend
