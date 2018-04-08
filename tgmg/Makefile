# Templated geometric multigrid library (TGMG)
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : December 6th 2015

#Load the common configuration file
include ../config.mk

#List of the common source files
tgmg_src=tgmg_config.hh tgmg.hh tgmg.cc tgmg_predict.hh
execs=poisson

#Makefile rules
all: $(execs)

depend: $(src)
	$(cxx) $(cflags) -MM $(src) >Makefile.dep

poisson: poisson.cc $(tgmg_src)
	$(cxx) $(cflags) -o $@ $<

%.o: %.cc
	$(cxx) $(cflags) -c $<

tgmg.cc: tgmg.in.cc tgmg_compile.pl
	perl tgmg_compile.pl

clean:
	rm -f $(execs)

.PHONY: all clean
