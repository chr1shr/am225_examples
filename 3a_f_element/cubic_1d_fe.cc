#include <cstdlib>

#include "cubic_1d_fe.hh"
#include "blas.h"

void cubic_1d_fe::init_const() {
    for(int i=0;i<3*n;i++) f[i]=1.;
    assemble_b();
}

void cubic_1d_fe::init_slope() {
    double xx=1;
    for(int i=0;i<3*n;i++,xx+=h) f[i]=1.5-xx;
    assemble_b();
}

void cubic_1d_fe::init_mms() {
    const double o=5*M_PI;
    double xx=1;
    for(int i=0;i<3*n;i++,xx+=h) f[i]=exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx));
    assemble_b();
}

void cubic_1d_fe::mul_A(double *in,double *out) {
    int i,j,k;
    const double B[16]={37/30.,-63/50.,9/20.,-13/120.,
                        -63/50.,18/5.,-99/40.,9/20.,
                        9/20.,-99/40.,18/5.,-63/50.,
                        -13/120.,9/20.,-63/50.,37/30.},
                 C[16]={17/30.,-51/80,3/8.,-13/80.,
                        -51/80.,27/8.,-297/80.,39/40.,
                        3/8.,-297/80.,297/40.,-327/80.,
                        -13/80.,39/40.,-327/80.,131/80.};
    for(i=0;i<3*n;i++) out[i]=0.;

    for(k=0;k<3*n;k+=3) for(i=(k==0?1:0);i<4;i++)
            for(j=0;j<4;j++) out[-1+k+i]+=((3*k+1./h)*B[i+4*j]+C[i+4*j])*in[-1+k+j];

    for(i=0;i<3*n;i++) out[i]*=1./h;
}

void cubic_1d_fe::assemble_b() {
    int i,j,k;
    const double D[16]={8/35.,99/560.,-9/140.,19/560.,
                        99/560.,81/70.,-81/560.,-9/140.,
                        -9/140.,-81/560.,81/70.,99/560.,
                        19/560.,-9/140.,99/560.,8/35.};

    for(i=0;i<3*n;i++) b[i]=0.;

    for(int k=0;k<3*n;k+=3) for(int i=(k==0?1:0);i<4;i++)
        for(int j=0;j<4;j++) b[-1+k+i]+=D[i+4*j]*f[k+j];

    for(i=0;i<3*n;i++) b[i]*=h;
}

void cubic_1d_fe::print(FILE *fp) {
    double xx=1;
    fputs("1 0\n",fp);
    for(int i=0;i<3*n;i++,xx+=h) fprintf(fp,"%g %g\n",xx,x[i]);
}

double cubic_1d_fe::l2_norm_mms() {
    double l2=0.,xx=1.+h;
    for(int i=0;i<3*n-1;i++) {
        xx+=h;
        l2+=mms_dsq(xx,x[i]);
    }
    l2+=0.5*mms_dsq(2.,x[3*n-1]);
    return sqrt(h*l2);
}

void cubic_1d_fe::print(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    print(fp);
    fclose(fp);
}
