# Load the common configuration file
include ../config.mk

iflags=-I../tgmg
lflags=-L.

objs=common.o fluid_2d.o mgs_fem.o
src=$(patsubst %.o,%.cc,$(objs))
execs=fluid_test

all:
	$(MAKE) -C ../tgmg
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libf2d.a: $(objs)
	rm -f libf2d.a
	ar rs libf2d.a $^

fluid_test: fluid_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) libf2d.a

clean-all:
	$(MAKE) -C ../tgmg clean
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean all executables depend
