# Include the configuration file from the parent directory
include ../config.mk

# Lists of executables
execs=cpp_example1 cpp_example2 cpp_example3 cpp_example4 cpp_example5 \
	  diffuse false_sharing happy openmp_example1 openmp_example2 \
	  openmp_example2b openmp_example3 poly_timing race_condition \
	  race_condition2 ridders_array sum_atomic sum_critical sum_reduction

all: $(execs)

# A Makefile target to remove all the built files
clean:
	rm -f $(execs)

# Makefile rule to build an executable. Note $@ is a shorthand for the file to
# be built and $^ is a shorthand for all the dependencies. (Can also use $< to
# mean the first dependency only.)
cpp_example1: cpp_example1.cc
	$(cxx) $(cflags) -o $@ $^

cpp_example2: cpp_example2.cc
	$(cxx) $(cflags) -o $@ $^

cpp_example3: cpp_example3.cc
	$(cxx) $(cflags) -o $@ $^

cpp_example4: cpp_example4.cc
	$(cxx) $(cflags) -o $@ $^

cpp_example5: cpp_example5.cc
	$(cxx) $(cflags) -o $@ $^

diffuse: diffuse.cc
	$(cxx) $(cflags) -o $@ $^

false_sharing: false_sharing.cc
	$(cxx) $(cflags) -o $@ $^

happy: happy.cc
	$(cxx) $(cflags) -o $@ $^

openmp_example1: openmp_example1.cc
	$(cxx) $(cflags) -o $@ $^

openmp_example2: openmp_example2.cc
	$(cxx) $(cflags) -o $@ $^

openmp_example2b: openmp_example2b.cc
	$(cxx) $(cflags) -o $@ $^

openmp_example3: openmp_example3.cc
	$(cxx) $(cflags) -o $@ $^

poly_timing: poly_timing.cc
	$(cxx) $(cflags) -o $@ $^

race_condition: race_condition.cc
	$(cxx) $(cflags) -o $@ $^

race_condition2: race_condition2.cc
	$(cxx) $(cflags) -o $@ $^

ridders_array: ridders_array.cc
	$(cxx) $(cflags) -o $@ $^

sum_atomic: sum_atomic.cc
	$(cxx) $(cflags) -o $@ $^

sum_critical: sum_critical.cc
	$(cxx) $(cflags) -o $@ $^

sum_reduction: sum_reduction.cc
	$(cxx) $(cflags) -o $@ $^

# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean
