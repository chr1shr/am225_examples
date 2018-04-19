# Load the common configuration file
include ../config.mk

iflags=-I../tgmg
lflags=-L.

objs=common.o levelset.o
src=$(patsubst %.o,%.cc,$(objs))
execs=ls_test

all:
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libf2d.a: $(objs)
	rm -f liblsm.a
	ar rs liblsm.a $^

ls_test: ls_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -llsm

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) liblsm.a

.PHONY: clean all executables depend
