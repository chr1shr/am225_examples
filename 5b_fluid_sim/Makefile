# 2D fluid simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : May 30th 2014

# Load the common configuration file
include ../../config.mk

iflags=-I../../tgmg -I../../shared
lflags=-L../../shared -L.

objs=common.o fluid_2d.o multisetup.o visco_impl.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mg_test conv_test fluid_test edges tconv

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../tgmg
	$(MAKE) executables
	
executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) -Imulti >Makefile.dep

include Makefile.dep

libf2d.a: $(objs)
	rm -f libf2d.a
	ar rs libf2d.a $^

fluid_test: fluid_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

tconv: tconv.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

conv_test: conv_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

mg_test: mg_test.cc multisetup.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^

edges: edges.cc common.o
	$(im_cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lgpmtx $(im_lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) libf2d.a

clean-all:
	$(MAKE) -C ../../shared clean
	$(MAKE) -C ../../tgmg clean
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean all executables depend
