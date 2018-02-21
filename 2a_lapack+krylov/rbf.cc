#include "rbf.hh"
#include "lp_solve.hh"
#include "file_output.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

/** Initializes the radial basis function class, setting constants
 * and allocating memory.
 * \param[in] n_ the number of points.
 * \param[in] type_ the type of radial function to use. */
rbf::rbf(int n_,int type_) : conj_grad(n_), n(n_), type(type_), lsq(1.), ilsq(1.),
    px(new double[n]), py(new double[n]), pf(new double[n]), rs(new double[n]),
    A(NULL), Apre(NULL), ipiv(NULL), An(NULL), phi_index(NULL), phi_val(NULL) {}

/** The class destructor frees the dynamically allocated memory. */
rbf::~rbf() {

    // Remove sparse matrix table if present
    if(An!=NULL) {
        delete [] phi_val;
        delete [] phi_index;
        delete [] An;
    }

    // Remove preconditioning table if present
    if(Apre!=NULL) {
        delete [] Apre;
        delete [] ipiv;
    }

    // Remove dense matrix memory if present
    if(A!=NULL) delete [] A;

    // Remove sample point information
    delete [] rs;
    delete [] pf;
    delete [] py;
    delete [] px;
}

/** Initializes the points to be randomly distributed in the square [-1,1]^2.
 */
void rbf::init_random() {

    // Create random positions and resort y positions
    for(int i=0;i<n;i++) {
        px[i]=urand();
        py[i]=urand();
    }
    std::sort(py,py+n);

    // Create function values
    for(int i=0;i<n;i++) {
        pf[i]=exp(-2*(px[i]*px[i]+py[i]*py[i]));
        rs[i]=pf[i];
    }
}

/** Computes the weights in the RBF interpolant using the direct LAPACK solver. */
void rbf::solve_weights_lapack() {
    assemble_matrix();
    solve_sym_matrix(n,A,rs);
}

/** Solve weights using conjugate gradient method.
 * \param[in] bls the block size to use in the block Jacobi preconditioner. If
 *                set to zero, then no preconditioning is used. */
void rbf::solve_weights_conj_grad(int bls,bool verbose) {
    make_table();copy(pf,b);
    if(bls>0) {
        preconditioning_table(bls);
        solve_pre(verbose);
    } else solve(verbose);
    copy(x,rs);
}

/** Sets up a table for use with preconditioning. This is a block Jacobi
 * preconditioner. The routine computes the diagonal blocks, and performs
 * the LU decomposition on each so that the inverse can be computed quickly
 * \param[in] bls_ the block size. */
void rbf::preconditioning_table(int bls_) {

    // Remove any previous preconditioning information
    if(Apre==NULL) {
        delete [] ipiv;delete [] Apre;
    }

    // Set constants and allocate memory
    bls=bls_;lbls=n%bls;
    int fb=n/bls,k;
    Apre=new double[fb*bls*bls+lbls*lbls];
    ipiv=new int[n];

    // Compute each diagonal block and perform the LU decomposition
    double *Ap=Apre;
    for(k=0;k<=n-bls;k+=bls,Ap+=bls*bls) {
        fill_matrix_entries(Ap,k,bls);
        factor_sym(bls,Ap,ipiv+k);
    }

    // Compute the final block, which may be of smaller size
    if(k<n) {
        fill_matrix_entries(Ap,k,lbls);
        factor_sym(lbls,Ap,ipiv+k);
    }
}

/** Assembles the RBF linear system in an array suitable for passing to LAPACK
 * functions. Due to symmetry, only the upper triangular part of the matrix is
 * filled in. */
void rbf::fill_matrix_entries(double *Ap,int k,int b) {
    double *plx=px+k,*ply=py+k;
    double dx,dy;
    for(int j=0;j<b;j++) {

        // Deal with off-diagonal elements, only filling in the upper
        // triangular part
        for(int i=0;i<j;i++) {
            dx=plx[i]-plx[j];
            dy=ply[i]-ply[j];
            Ap[i+b*j]=phi(dx*dx+dy*dy);
        }

        // Set diagonal element
        Ap[(b+1)*j]=phi(0);
    }
}

/** Outputs the points, their function values, and their weights in the RBF
 * interpolant.
 * \param[in] filename the name of the file to write to. */
void rbf::output_points(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    for(int i=0;i<n;i++) fprintf(fp,"%g %g %g %g\n",px[i],py[i],pf[i],rs[i]);
    fclose(fp);
}

/** Outputs the radial basis function interpolant.
 * \param[in] filename the name of the file to write to.
 * \param[in] q the number of grid points in each direction.
 * \param[in] L the maximum coordinate extent. */
void rbf::output_interpolant(const char* filename,int q,double L) {
    const double h=2.*L/(q-1);
    double *fld=new double[q*q],s,delx,dely,x,y;

    // Loop over the interpolation grid
    for(int j=0;j<q;j++) {
        y=-L+h*j;
        for(int i=0;i<q;i++) {
            x=-L+h*i;

            // Assemble RBF interpolant
            s=0;
            for(int k=0;k<n;k++) {
                delx=x-px[k];
                dely=y-py[k];
                s+=rs[k]*phi(delx*delx+dely*dely);
            }
            fld[i+q*j]=s;
        }
    }

    // Call library routine to output Gnuplot field
    gnuplot_output(filename,fld,q,q,-L,L,-L,L);
    delete [] fld;
}

/** Calculates the radial function.
 * \param[in] rsq the squared radial coordinate. */
double rbf::phi(double rsq) {
    if(type==0) {

        // Gaussian
        return exp(-rsq*ilsq);
    } else if(type==1) {

        // Polyharmonic spline
        return rsq==0?0:0.5*rsq*log(rsq);
    } else {

        // Wendland compactly supported function
        // (1-r)_+^8 (32r^3+25r^2+8r+1)
        if(rsq>lsq) return 0;
        rsq*=ilsq;
        double r=sqrt(rsq),mr=1-r;
        mr*=mr;mr*=mr;
        return mr*mr*(1+8*r+rsq*(25+32*r));
    }
}

/** Computes a uniformly distributed number between -1 and 1. */
double rbf::urand() {
    return -1.+(2./RAND_MAX)*double(rand());
}

/** Makes a table of non-zero entries in the RBF linear system, for the case of
 * a compactly supported radial function. */
void rbf::make_table() {
    if(An!=NULL) return;
    int c=0;
    double rsq,dx,dy;
    An=new int[n];

    // Count the number of non-zero terms
    for(int j=0;j<n;j++) {
        An[j]=0;
        for(int i=0;i<j;i++) {
            dx=px[i]-px[j];
            dy=py[i]-py[j];
            rsq=dx*dx+dy*dy;
            if(rsq<lsq) An[j]++;
        }
        c+=An[j];
    }

    // Allocate and store the sparse matrix table
    phi_index=new int[c];
    phi_val=new double[c];
    int *p_index=phi_index;
    double *p_val=phi_val;
    for(int j=0;j<n;j++) {
        for(int i=0;i<j;i++) {
            dx=px[i]-px[j];
            dy=py[i]-py[j];
            rsq=dx*dx+dy*dy;
            if(rsq<lsq) {
                *(p_index++)=i;
                *(p_val++)=phi(rsq);
            }
        }
    }
}

/** Performs black-box multiplication. */
void rbf::mul_A(double *in,double *out) {

    // Deal with diagonal terms
    for(int j=0;j<n;j++) out[j]=in[j];

    // Deal with off-diagonal terms using pre-computed table
    int i,j,*p_index=phi_index,*p2=p_index;
    double *p_val=phi_val;
    for(j=0;j<n;j++) {
        p2=p_index+An[j];
        while(p_index<p2) {
            i=*(p_index++);
            out[j]+=*p_val*in[i];
            out[i]+=*(p_val++)*in[j];
        }
    }
}

/** Performs preconditioning with the block Jacobi preconditioner M.
 * \param[in] in the input vector.
 * \param[in] out the output vector, after applying the M^{-1} matrix. */
void rbf::M_inv(double *in,double *out) {
    copy(in,out);
    int k;
    double *Ap=Apre;

    // Invert each block using the previously computed LU decompositions
    for(k=0;k<=n-bls;k+=bls,Ap+=bls*bls) factor_sym_solve(bls,Ap,ipiv+k,out+k);
    if(k<n) factor_sym_solve(lbls,Ap,ipiv+k,out+k);
}

/** Computes and prints the eigenvalues. It also prints the condition number */
void rbf::eigenvalues() {

    // Compute eigenvalues
    assemble_matrix();
    double *evals=new double[n];
    evals_sym_matrix(n,A,evals);
    double kappa=evals[n-1]/(*evals);

    // Compute condition number
    double emin=fabs(*evals),emax=emin,e;
    for(double *ep=evals+1;ep<evals+n;ep++) {
        e=fabs(*ep);
        if(e>emax) emax=e;
        if(e<emin) emin=e;
    }
    printf("# Condition number %g\n",emax/emin);

    // Print the eigenvalues
    for(int i=0;i<n;i++) printf("%g\n",evals[i]);
    delete [] evals;
}
