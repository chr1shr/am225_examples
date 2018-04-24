# Load the common configuration file
include ../config.mk

iflags=-I../tgmg
lflags=-L.

objs=common.o levelset.o fmm.o
src=$(patsubst %.o,%.cc,$(objs))
execs=ls_test fmm_test

all:
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

liblsm.a: $(objs)
	rm -f liblsm.a
	ar rs liblsm.a $^

ls_test: ls_test.cc liblsm.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -llsm

fmm_test: fmm_test.cc liblsm.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -llsm

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) liblsm.a

.PHONY: clean all executables depend
