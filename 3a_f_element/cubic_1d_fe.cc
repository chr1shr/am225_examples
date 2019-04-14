#include <cstdlib>

#include "cubic_1d_fe.hh"
#include "blas.h"

/** Initializes the source function to be a constant. */
void cubic_1d_fe::init_const() {
    for(int i=0;i<=3*n;i++) f[i]=1.;
    assemble_b();
}

/** Initializes the source function to be a linear slope, f(x)=1.5-x. */
void cubic_1d_fe::init_slope() {
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=xx-1.5;
    assemble_b();
}

/** Initializes the source function so that the solution will match
 * a manufactured solution, u(x)=exp(1-x)*sin(5*pi*x). */
void cubic_1d_fe::init_mms() {
    const double o=5*M_PI;
    double xx=1;
    for(int i=0;i<=3*n;i++,xx+=h) f[i]=-exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx));
    assemble_b();
}

/** Prints the stiffness matrix as a text array. */
void cubic_1d_fe::print_matrix() {
    int i,j;

    // Allocate zero vector, and workspace to compute a matrix product
    double *r=new double[6*n],*s=r+3*n;
    for(i=0;i<3*n;i++) r[i]=0.;

    for(int i=0;i<3*n;i++) {

        // Apply the black box matrix multiplication routine to a unit vector,
        // in order to extract a column of matrix entries
        r[i]=1.;mul_A(r,s);r[i]=0.;

        // Print a row of matrix entries. This assumes the matrix is symmetric
        // (as required for conjugate gradient) so that the row<->column switch
        // is permissible.
        for(j=0;j<3*n-1;j++) printf("%g ",s[j]);
        printf("%g\n",s[3*n-1]);
    }
    delete [] r;
}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_fe::mul_A(double *in,double *out) {
    int i,j,k;

    // Pre-computed integrals of derivatives of Lagrange polynomials, which are
    // required to construct the stiffness matrix
    const double B[16]={37/30.,-63/40.,9/20.,-13/120.,
                        -63/40.,18/5.,-99/40.,9/20.,
                        9/20.,-99/40.,18/5.,-63/40.,
                        -13/120.,9/20.,-63/40.,37/30.},
                 C[16]={17/40.,-51/80.,3/8.,-13/80.,
                        -51/80.,27/8.,-297/80.,39/40.,
                        3/8.,-297/80.,297/40.,-327/80.,
                        -13/80.,39/40.,-327/80.,131/40.};

    // Set the output vector to initially be zero
    for(i=0;i<3*n;i++) out[i]=0.;

    // Loop over each interval, and compute the contribution from
    // each
    for(k=0;k<3*n;k+=3) for(i=(k==0?1:0);i<4;i++)
        for(j=(k==0?1:0);j<4;j++)
            out[-1+k+i]+=((k+1./h)*B[i+4*j]+C[i+4*j])*in[-1+k+j];
}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void cubic_1d_fe::assemble_b() {
    int i,j,k;

    // Pre-computed integrals of Lagrange polynomial products
    const double D[16]={8/35.,99/560.,-9/140.,19/560.,
                        99/560.,81/70.,-81/560.,-9/140.,
                        -9/140.,-81/560.,81/70.,99/560.,
                        19/560.,-9/140.,99/560.,8/35.};

    // Clear the source function
    for(i=0;i<3*n;i++) b[i]=0.;

    // Loop over each interval, and compute the contributions
    // from each Lagrange polynomial pair
    for(k=0;k<3*n;k+=3) for(i=(k==0?1:0);i<4;i++)
        for(j=0;j<4;j++) b[-1+k+i]+=D[i+4*j]*f[k+j];

    // Normalize the results, and add in the Neumann condition to the last
    // entry
    for(i=0;i<3*n;i++) b[i]*=h;
    b[3*n-1]+=2*g;
}

/** Prints the solution.
 * \param[in] fp a file handle to write to. */
void cubic_1d_fe::print(FILE *fp) {
    double xx=1+h;
    fprintf(fp,"1 0 0 %g\n",*f);
    for(int i=0;i<3*n;i++,xx+=h) fprintf(fp,"%g %g %g %g\n",xx,x[i],b[i],f[i+1]);
}

/** Prints the solution.
 * \param[in] filename the name of the file to write to. */
void cubic_1d_fe::print(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    print(fp);
    fclose(fp);
}

/** Computes the L2 norm between the numerical solution and the true
 * manufactured solution. The routine uses the trapezoid rule to evaluate the
 * integral.
 * \return The L2 norm. */
double cubic_1d_fe::l2_norm_mms() {
    double l2=0.,xx=1.;
    for(int i=0;i<3*n-1;i++) {
        xx+=h;
        l2+=mms_dsq(xx,x[i]);
    }

    // Add the contribution at the last point, including a factor of 1/2 due to
    // the trapezoid rule. Note that there is no contribution at x=1, since the
    // numerical solution is zero there, and hence matches the manufactured
    // solution perfectly.
    l2+=0.5*mms_dsq(2.,x[3*n-1]);
    return sqrt(h*l2);
}
