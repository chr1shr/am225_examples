# Include the configuration file
include ../config.mk

# Flags for including the RK4 solver
ode_iflags=-I../1a_ode_solvers
cm_objs=../1a_ode_solvers/sol_base.o ../1a_ode_solvers/sol_rk4.o

# Lists of files to be built
objs=v_visco.o traffic.o car_model.o
src=$(patsubst %.o,%.cc,$(objs))
execs=t_solve cm_solve vv_solve

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

t_solve: t_solve.cc traffic.o
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

vv_solve: vv_solve.cc v_visco.o
	$(cxx) $(cflags) -o $@ $^ $(ode_iflags) $(lp_lflags)

cm_solve: cm_solve.cc car_model.o $(cm_objs)
	$(cxx) $(cflags) -o $@ $^ $(ode_iflags) $(lp_lflags)

.PHONY: clean depend executables
