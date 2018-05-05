#include <cstdio>
#include <cmath>
#include <cstring>

#include "sol_euler.hh"
#include "sol_improv_e.hh"
#include "sol_imp_onestep.hh"
#include "sol_sep.hh"
#include "orbit.hh"

int main(int argc,char **argv) {

    // Number of timesteps
    const int steps=1440;

    // The semi-major axis of the orbit
    const double a=1.;

    // The eccentricity of the orbit
    const double e=1/3.;

    // The simulation duration, currently set to three complete orbits
    const double du=6*M_PI*a*sqrt(a);

    // Check for one command-line argument
    if(argc!=2) {
        fputs("Syntax: ./orb_solve <method_type>\n\n"
              "Solves an elliptical orbit problem with a variety of non-symplectic (NS) and\n"
              "symplectic (S) methods.\n\n"
              "Types:\n"
              "'fwd' - forward Euler, NS\n"
              "'bck' - backward Euler, NS\n"
              "'ime' - improved Euler, NS\n"
              "'mpt' - implicit midpoint method, S\n"
              "'fos' - first-order symplectic method, S\n"
              "'ru3' - Ruth's third-order method, S\n"
              "'ru4' - Ruth's fourth-order method, S\n",stderr);
        return 1;
    }

    // Check the method type, and call the corresponding integration routine
    if(strcmp(argv[1],"fwd")==0) {

        // Explicit forward Euler
        orb_euler o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"bck")==0) {

        // Backward Euler
        orb_back o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"ime")==0) {

        // Improved Euler
        orb_improv_e o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"mpt")==0) {

        // Implicit midpoint
        orb_imp_mid o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"fos")==0) {

        // First-order symplectic
        orb_fo_symplectic o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"ru3")==0) {

        // Ruth's third-order method
        orb_ruth3 o(a,e);o.solve_fixed(du,steps,true);
    } else if(strcmp(argv[1],"ru4")==0) {

        // Ruth's fourth-order method
        orb_ruth4 o(a,e);o.solve_fixed(du,steps,true);
    } else fputs("Invalid method type\n",stderr);
}
