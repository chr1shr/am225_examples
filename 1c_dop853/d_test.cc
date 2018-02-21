#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "dop853.hh"
#include "functions.hh"

/** The number of components in the ODE system. */
const int ns=2;

/** The number of dense outputs to store. */
const int n_do=2001;

/** The initial integration time. */
const double tstart=0;

/** The final integration time. */
const double tend=20;

int main() {
	double wmem[n_do*ns],*ws[n_do],*wp,**wsp=ws,t;

	// Initialize pointers for saving the dense output interpolations
	for(wp=wmem;wp<wmem+n_do*ns;wp+=ns) *(wsp++)=wp;

	// Construct the test class, set the level of accuracy, and carry out
	// the ODE solution
	brusselator o;
	o.atoli=1e-6;
	o.solve(tstart,tend,n_do,ws);

	// Output the saves
	FILE *fp=safe_fopen("dt.out2");
	const double dt=(tend-tstart)/(n_do-1);
	for(int i=0;i<n_do;i++) {
		t=i*dt;
		fprintf(fp,"%.15g %.15g %.15g\n",t,ws[i][0],
			ws[i][1]);
//		fprintf(fp,"%.15g %.15g %.15g %.15g %.15g\n",t,ws[i][0],
//			ws[i][1],ws[i][0]-o.sol0(t),ws[i][1]-o.sol1(t));
	}
	fclose(fp);
}
