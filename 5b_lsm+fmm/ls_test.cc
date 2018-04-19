#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "levelset.hh"

const char fn[]="lsm.out";

int main() {

    // Create the output directory for storing the simulation frames
    mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Construct the simulation class, setting the number of gridpoints, the
    // periodicity, and physical constants
    levelset ls(256,256,true,true,-1,1,-1,1,fn);

    // Initialize the tracers, and set the timestep based on multiplying the
    // maximum allowable by a padding factor
    ls.initialize(0,0.2);

    // Run the simulation for a specified duration, outputting snapshots at
    // regular intervals
    ls.solve(2*M_PI,180);
}
