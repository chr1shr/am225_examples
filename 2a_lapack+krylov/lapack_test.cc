#include <cstdio>
#include <cstdlib>

#include "lp_solve.hh"

int main() {
	double A[4]={1,1,2,1};
	double b[2]={3,7};

	// Print the matrix and the right hand side
	printf("A=[%6g %6g ]\n  [%6g %6g ]\n\n",*A,A[2],A[1],A[3]);
	printf("b=[%6g ]\n  [%6g ]\n\n",*b,b[1]);

	// Call routine to solve the matrix
	solve_matrix(2,A,b);

	// Print the solution
	printf("x=[%6g ]\n  [%6g ]\n",*b,b[1]);
}
