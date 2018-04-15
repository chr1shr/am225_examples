#include <cstdio>
#include <cstring>

#include "ho_transport.hh"

int main(int argc,char **argv) {

    // Print a syntax message if a command-line argument isn't provided
    if(argc!=2) {
        fputs("Syntax: ./ho_solve <type>\n\n"
              "<type> can be: 'gd' for the Godunov method\n"
              "               'lw' for the Lax-Wendroff method\n"
              "               'bw' for the Beam-Warming method\n"
              "               'mm' for the minmod limiter method\n"
              "               'sb' for the superbee limiter method\n"
              "               'e2' for the ENO2 method\n",stderr);
        return 1;
    }

    // Read the integration type
    int type;
    if(strcmp(argv[1],"gd")==0) type=0;
    else if(strcmp(argv[1],"lw")==0) type=1;
    else if(strcmp(argv[1],"bw")==0) type=2;
    else if(strcmp(argv[1],"mm")==0) type=3;
    else if(strcmp(argv[1],"sb")==0) type=4;
    else if(strcmp(argv[1],"e2")==0) type=5;
    else {
        fputs("Unknown integration type\n",stderr);
        return 1;
    }

    // Number of snapshots to output
    const int snaps=20;

    // Number of gridpoints
    const int m=256;

    // Integration timestep safety factor
    const double sf=0.2;

    // Create the transport simulation, and set up the initial condition
    ho_transport ho(m,1.0);
//    ho.init_exp_sine();
    ho.init_step_function();

    // Integrate and save the solution snapshots to file
    char buf[32];
    sprintf(buf,"ho_%s.out",argv[1]);
    ho.solve(buf,snaps,1.0,sf,type);
}
