#include "lp_solve.hh"

#include <cstdio>
#include <cstdlib>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    void dgesv_(int *n,int *nrhs,double *a,int *lda,int *ipivot,
            double *b,int *ldb,int *info);
    void dgbsv_(int *n,int *kl,int *ku,int *nrhs,double *ab,int *ldab,int *ipiv,
            double *b,int *ldb,int *info);
    void dsysv_(char* uplo,int *n,int *nrhs,double *a,int *lda,int *ipivot,
            double *b,int *ldb,double *work,int *lwork,int *info);
    void dsyev_(char* jobz,char* uplo,int *n,double *a,int *lda,double *w,
            double *work,int *lwork,int *info);
    void dsytrf_(char* uplo,int *n,double *a,int *lda,int *ipivot,double *work,
            int *lwork,int *info);
    void dsytrs_(char* uplo,int *n,int *nrhs,double *a,int *lda,int *ipivot,
            double *b,int *ldb,int *info);
}

/** Prints a generic error message in cases when LAPACK reports an error. */
void lapack_fail() {
    fputs("LAPACK routine failed\n",stderr);
    exit(1);
}

/** Solves the matrix system Ax=b.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void solve_matrix(int n,double *A,double *x) {

    // Create the temporary memory that LAPACK needs
    int info,nrhs=1,*ipiv=new int[n];

    // Make call to LAPACK routine
    dgesv_(&n,&nrhs,A,&n,ipiv,x,&n,&info);

    // Remove temporary memory
    delete [] ipiv;

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0) lapack_fail();
}

/** Solves the matrix system Ax=b for the case when A is symmetric.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void solve_sym_matrix(int n,double *A,double *x) {

    // Create the temporary memory that LAPACK needs
    char uplo='u';
    int info,nrhs=1,*ipiv=new int[n],lwork=64*n;
    double *work=new double[lwork];

    // Make LAPACK call
    dsysv_(&uplo,&n,&nrhs,A,&n,ipiv,x,&n,work,&lwork,&info);

    // Remove temporary memory
    delete [] work;
    delete [] ipiv;

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0) lapack_fail();
}

/** Solves the matrix system Ax=b for the case when A is banded.
 * \param[in] n the dimension of the matrix.
 * \param[in] (kl,lu) the number of subdiagonals and superdiagonals.
 * \param[in] A the matrix terms; see DGBSV documentation for full details.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void solve_banded_matrix(int n,int kl,int ku,double *A,double *x) {

    // Create the temporary memory that LAPACK needs
    int info,nrhs=1,*ipiv=new int[n],ldab=2*kl+ku+1;

    // Make call to LAPACK routine
    dgbsv_(&n,&kl,&ku,&nrhs,A,&ldab,ipiv,x,&n,&info);
    if(info!=0) {
        fputs("LAPACK routine failed\n",stderr);
        exit(1);
    }

    // Remove temporary memory
    delete [] ipiv;
}

/** Finds the eigenvalues of a symmetric matrix A
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored.
 * \param[in] evals the eigenvalues. */
void evals_sym_matrix(int n,double *A,double *evals) {
    char jobz='n',uplo='u';
    int lwork=34*n,info;
    double *work=new double[lwork];

    // Make LAPACK call and remove temporary memory
    dsyev_(&jobz,&uplo,&n,A,&n,evals,work,&lwork,&info);
    delete [] work;

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0) lapack_fail();
}

/* Factors a symmetric matrix.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored. On exit the terms will contain the factorized
 *              matrix.
 * \param[in] ipiv the pivoting information of the factorization. */
void factor_sym(int n,double *A,int *ipiv) {

    // Create the temporary memory that LAPACK needs
    char uplo='u';
    int info,lwork=64*n;
    double *work=new double[lwork];

    // Perform LU decomposition
    dsytrf_(&uplo,&n,A,&n,ipiv,work,&lwork,&info);
    if(info!=0) lapack_fail();

    // Remove temporary memory
    delete [] work;
}

/** Solves a symmetric linear system using a previously factored matrix.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the factored matrix terms.
 * \param[in] ipiv the pivoting information.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void factor_sym_solve(int n,double *A,int *ipiv,double *x) {
    char uplo='u';
    int info,nrhs=1;
    dsytrs_(&uplo,&n,&nrhs,A,&n,ipiv,x,&n,&info);

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0) lapack_fail();
}
