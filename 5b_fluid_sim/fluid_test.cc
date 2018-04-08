#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "fluid_2d.hh"
#include "common.hh"

const char fn[]="ftest";

int main() {
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	fluid_2d f2d(128,128,true,true,-1,1,-1,1,500,1,25.,128,true,fn);
	f2d.init_fields();
	f2d.solve(0,1,200);
}
