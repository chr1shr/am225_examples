#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "fluid_2d.hh"
#include "common.hh"

const char fn[]="ftest";

int main() {
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    unsigned int fflags=7;
	fluid_2d f2d(256,256,true,true,-1,1,-1,1,0.002,1.,fflags,fn);
	f2d.initialize(512,0.6);
	f2d.solve(10,200);
}
