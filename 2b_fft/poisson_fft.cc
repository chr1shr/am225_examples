#include <cstdio>
#include <cmath>

#include "file_output.hh"
#include "poisson_fft.hh"

poisson_fft::poisson_fft(int n_)
    : n(n_), nn(n*n), h(2./(n+1)),
    f(fftw_alloc_real(nn)), v(fftw_alloc_real(nn)), w(fftw_alloc_real(nn)),
    plan_fwd(fftw_plan_r2r_2d(n,n,f,w,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_bck(fftw_plan_r2r_2d(n,n,w,v,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)) {}

poisson_fft::~poisson_fft() {
    fftw_destroy_plan(plan_bck);
    fftw_destroy_plan(plan_fwd);
    fftw_free(w);
    fftw_free(v);
    fftw_free(f);
}

void poisson_fft::init() {
    for(int j=0;j<n;j++) {
        double y=-1+(j+1)*h;
        for(int i=0;i<n;i++) {
            double x=-1+(i+1)*h;
            f[i+n*j]=fabs(x)<0.5&&fabs(y)<0.5?(x>0?-1:1):0;
        }
    }
}

void poisson_fft::solve() {
    const double nor=1./(2*(n+1));
    fftw_execute(plan_fwd);

    for(int j=0;j<n;j++) {
        double lamj=lambda(j);
        for(int i=0;i<n;i++) {
            w[i+n*j]*=h*h*nor*nor/(lambda(i)+lamj);
        }
    }
    fftw_execute(plan_bck);
}

double poisson_fft::lambda(int j) {
    return 2*(1-cos(M_PI/(n+1)*(j+1)));
}

void poisson_fft::print() {
    for(int j=0;j<n;j++) {
        printf("%g",v[j]);
        for(int i=1;i<n;i++) {
            printf(" %g",v[j+i*n]);
        }
        putchar('\n');
    }
}

void poisson_fft::output(const char* filename,bool solution) {

    double *fld=new double[(n+2)*(n+2)],*f2=fld+(n+3);
    for(int j=0;j<(n+2)*(n+2);j++) fld[j]=0;

    for(int j=0;j<n;j++) for(int i=0;i<n;i++)
        f2[i+j*(n+2)]=solution?v[i+j*n]:f[i+j*n];

    gnuplot_output(filename,fld,n+2,n+2,-1,1,-1,1);

    delete [] fld;
}

void poisson_fft::init_mms() {
    for(int j=0;j<n;j++) {
        double y=-1+(j+1)*h;
        for(int i=0;i<n;i++) {
            double x=-1+(i+1)*h;
            f[i+n*j]=exp(x)*(3-y*y-4*x*(y*y-1)-x*x*(1+y*y));
        }
    }
}

double poisson_fft::l2_error_mms() {
    double l2=0,del;
    for(int j=0;j<n;j++) {
        double y=-1+(j+1)*h;
        for(int i=0;i<n;i++) {
            double x=-1+(i+1)*h;
            del=v[i+n*j]-exp(x)*(1-x*x)*(1-y*y);
            l2=del*del;
        }
    }
    return sqrt(h*h*l2);
}
