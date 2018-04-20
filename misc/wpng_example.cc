#include "write_png.hh"

#include <cmath>

// Function for making a rainbow palette, identical to that defined in HW2
// question 3.
double f(double t) {return 0.45*(1+cos(t));}

int main() {

    // Size of ouptut image
    const int m=512,n=512,mn=m*n;

    // Grid spacing
    const double h=2./n;

    // Allocate memory for the three color channels per gridpoint
    double *z=new double[3*mn],*zp=z,fr,t,x,y;

    // Loop over the pixels and set the colors according to a radial pattern
    for(int j=0;j<n;j++) {
        y=-1+h*(j+0.5);
        for(int i=0;i<m;i++) {
            x=-1+h*(i+0.5);
            t=atan2(y,x);
            fr=sin(6*M_PI*sqrt(x*x+y*y));
            *(zp++)=f(t+fr);
            *(zp++)=f(t+fr-2*M_PI/3.);
            *(zp++)=f(t+fr+2*M_PI/3.);
        }
    }

    // Call routine to write PNG file, and free the dynamically allocated
    // memory
    write_png("test.png",m,n,z,0.,1.);
    delete [] z;
}
