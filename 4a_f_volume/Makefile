# Include the configuration file
include ../config.mk

# Flags for including the RK4 solver
ode_iflags=-I../1a_ode_solvers
cm_objs=../1a_ode_solvers/sol_base.o ../1a_ode_solvers/sol_rk4.o car_model.o

# Lists of files to be built
objs=diffuse.o ho_transport.o traffic.o
src=$(patsubst %.o,%.cc,$(objs)) car_model.cc
execs=d_solve d_cons_test ho_solve t_solve cm_solve

all:
	$(MAKE) -C ../1a_ode_solvers sol_base.o sol_rk4.o
	$(MAKE) executables

executables: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) $(ode_iflags) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $< $(ode_iflags)

ho_solve: ho_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

d_solve: d_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

d_cons_test: d_cons_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

t_solve: t_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

cm_solve: cm_solve.cc $(objs) $(cm_objs)
	$(cxx) $(cflags) -o $@ $^ $(ode_iflags) $(lp_lflags)

.PHONY: clean depend executables
