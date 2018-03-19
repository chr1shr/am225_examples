#include <cstdlib>

#include "cubic_1d_fe.hh"
#include "blas.h"

void cubic_1d_fe::init_const() {
    for(int i=0;i<=3*n;i++) f[i]=1.;
    assemble_b();
}

void cubic_1d_fe::init_slope() {
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=1.5-xx;
    assemble_b();
}

void cubic_1d_fe::init_mms() {
    const double o=5*M_PI;
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx));
    assemble_b();
}

void cubic_1d_fe::print_matrix() {
    int i,j;
    double *r=new double[6*n],*s=r+3*n;

    for(i=0;i<3*n;i++) r[i]=0.;

    for(int i=0;i<3*n;i++) {
        r[i]=1.;
        mul_A(r,s);
        r[i]=0.;
        for(j=0;j<3*n-1;j++) printf("%g ",s[j]);
        printf("%g\n",s[3*n-1]);
    }

    delete [] r;
}

void cubic_1d_fe::mul_A(double *in,double *out) {
    int i,j,k;
    const double B[16]={37/30.,-63/40.,9/20.,-13/120.,
                        -63/40.,18/5.,-99/40.,9/20.,
                        9/20.,-99/40.,18/5.,-63/40.,
                        -13/120.,9/20.,-63/40.,37/30.},
                 C[16]={17/40.,-51/80.,3/8.,-13/80.,
                        -51/80.,27/8.,-297/80.,39/40.,
                        3/8.,-297/80.,297/40.,-327/80.,
                        -13/80.,39/40.,-327/80.,131/40.};
    for(i=0;i<3*n;i++) out[i]=0.;

    for(k=0;k<3*n;k+=3) for(i=(k==0?1:0);i<4;i++)
            for(j=(k==0?1:0);j<4;j++) out[-1+k+i]+=((k+1./h)*B[i+4*j]+C[i+4*j])*in[-1+k+j];
}

void cubic_1d_fe::assemble_b() {
    int i,j,k;
    const double D[16]={8/35.,99/560.,-9/140.,19/560.,
                        99/560.,81/70.,-81/560.,-9/140.,
                        -9/140.,-81/560.,81/70.,99/560.,
                        19/560.,-9/140.,99/560.,8/35.};

    for(i=0;i<3*n;i++) b[i]=0.;

    for(k=0;k<3*n;k+=3) for(i=(k==0?1:0);i<4;i++)
        for(j=0;j<4;j++) b[-1+k+i]-=D[i+4*j]*f[k+j];

    for(i=0;i<3*n;i++) b[i]*=h;
    b[3*n-1]+=2*g;
}

void cubic_1d_fe::print(FILE *fp) {
    double xx=1+h;
    fprintf(fp,"1 0 0 %g\n",*f);
    for(int i=0;i<3*n;i++,xx+=h) fprintf(fp,"%g %g %g %g\n",xx,x[i],b[i],f[i+1]);
}

double cubic_1d_fe::l2_norm_mms() {
    double l2=0.,xx=1.;
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
