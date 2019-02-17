# Specify compiler and flags
cxx=g++
cflags=-Wall

# Rule to build object file
file_output.o: file_output.cc file_output.hh
	$(cxx) $(cflags) -c file_output.cc

# Rule to build executable
gp_test: gp_test.cc file_output.o
	$(cxx) $(cflags)  -o gp_test gp_test.cc file_output.o
