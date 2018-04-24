#include "fmm.hh"

int main() {

    fmm f(128,128,0,1,0,1);

    f.init_fields();

    f.fast_march();

    f.output("phi.out",0);
    f.output("ind.out",1);
}
