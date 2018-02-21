#ifndef LP_SOLVE_HH
#define LP_SOLVE_HH

void lapack_fail();
void solve_matrix(int n,double *A,double *x);
void solve_sym_matrix(int n,double *A,double *x);
void evals_sym_matrix(int n,double *A,double *eval);
void factor_sym(int n,double *A,int *ipiv);
void factor_sym_solve(int n,double *A,int *ipiv,double *x);

#endif
