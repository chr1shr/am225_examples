#include <cmath>

#include "file_output.hh"

int main() {

    // Physical dimensions of the grid
    const double ax=-1,bx=1,ay=-1,by=1;

    // Size of grid
    const int m=129,n=129;

    // Grid spacings
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

    // Allocate memory and fill the grid with a test function
    double *fld=new double[m*n],*fp=fld;
    for(int j=0;j<n;j++) {
        double y=ay+j*dy;
        for(int i=0;i<m;i++) {
            double x=ax+i*dx;
            *(fp++)=sin(M_PI*x)*cos(M_PI*(x+3*y));
        }
    }

    // Call output function and clean up memory
    gnuplot_output("test_out.gnu",fld,m,n,ax,bx,ay,by);
    delete [] fld;
}
