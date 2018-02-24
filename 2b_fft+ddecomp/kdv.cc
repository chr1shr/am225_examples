#include "kdv.hh"

#include <cstdlib>
#include <cmath>

void kdv::ff(double t_,double *in,double *out) {
    double re,im;
    const double ninv=1./n;
    for(int i=0;i<n;i++) der[i]=in[i];

    // Compute the first and third derivatives
    fftw_execute(plan1);
    for(int i=0;i<fftn;i++) {
        re=c1[i][0];
        im=c1[i][1];

	// First derivative Fourier spectrum
        c1[i][0]=-ninv*i*im;
        c1[i][1]=ninv*i*re;

	// Third derivative Fourier spectrum
        c2[i][0]=ninv*i*i*i*im;
        c2[i][1]=-ninv*i*i*i*re;
    }
    fftw_execute(plan2);
    fftw_execute(plan3);

    for(int i=0;i<n;i++) out[i]=-in[i]*der[i]-0.01*der3[i];
}

void kdv::init() {
    const double h=2*M_PI/n;
    double x;

    for(int i=0;i<n;i++) {
        x=i*h;
        q[i]=exp((cos(x)-1)*4)+0.5*exp((cos(x-2)-1)*6);
    }
}
