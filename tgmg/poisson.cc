// This is an example file for testing the multigrid code.

#include "tgmg.cc"

// Multisetup structure for Poisson problem
struct multisetup1 {
	/** Grid dimensions. */
	const int m;
	const int n;
	/** Total number of gridpoints. */
	const int mn;
	/** Periodicity in the x and y directions. */
	const bool x_prd;
	const bool y_prd;
	/** The mode to use for the Gauss-Seidel smoothing. (0=default) */
	const char gs_mode;
	/** Lower and upper limits in the x direction. */
	const double ax,bx;
	/** Lower and upper limits in the y direction. */
	const double ay,by;
	/** Grid spacings in the x and y directions. */
	const double dx,dy;
	/** Stencil entries. */
	const double fm,fm_inv,fex,fey,fc;
	/** Threshold on L_2 norm of residual to terminate the multigrid solve. */
	const double acc;
	/** A pointer to the solution vector. */
	double* const z;
	multisetup1(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_,double* const z_)
		: m(m_), n(n_), mn(m_*n_), x_prd(false), y_prd(false),
		gs_mode(1), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), fm(-4/(dx*dx)),
		fm_inv(-0.25*dx*dx), fex(1/(dx*dx)), fey(1/(dx*dx)), fc(0),
		acc(tgmg_accuracy(fm,1e4)), z(z_) {}
	/** Function to determine whether a grid point is on the edge or not.
	 */
	inline bool edge(int i,int ij) {return i==0||i==m-1||ij>=mn-m||ij<m;}
	/** Functions to specify the corner stencil entries. */
	inline double a_dl(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_dr(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_ul(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_ur(int i,int ij) {return edge(i,ij)?0:fc;}
	/** Functions to specify the vertical stencil entries. */
	inline double a_dc(int i,int ij) {return edge(i,ij)?0:fey;}
	inline double a_uc(int i,int ij) {return edge(i,ij)?0:fey;}
	/** Functions to specify the horizontal stencil entries. */
	inline double a_cl(int i,int ij) {return edge(i,ij)?0:fex;}
	inline double a_cr(int i,int ij) {return edge(i,ij)?0:fex;}
	/** Function to specify the central stencil entry (on the diagonal of
	 * the linear system). */
	inline double a_cc(int i,int ij) {return fm;}
	/** Function to multiply by the reciprocal of the central stencil
	 * entry. This is specified as a separate function for computational
	 * efficiency. */
	inline double inv_cc(int i,int ij,double v) {return fm_inv*v;}
	/** Calculates the ith component of the multiplication (A-D)z, needed
	 * in the Gauss--Seidel smoothing iteration. */
	inline double mul_a(int i,int ij) {
		return edge(i,ij)?0:fex*(z[ij+1]+z[ij-1])+fey*(z[ij+m]+z[ij-m]);
	}
};

int main() {
	const int m=1025,n=1025,mn=m*n;
	const double ax=-8,bx=8,ay=ax,by=bx;
	int i,j,ij;
	double *b=new double[mn],*z=new double[mn],x,y;
	multisetup1 msu(m,n,ax,bx,ay,by,z);
	tgmg<multisetup1,double,double> mg(msu,b,z);
	mg.verbose=3;

	tgmg_predict tp;

	// Set up the multigrid hierarchy
	mg.setup();

	for(int k=0;k<50;k++) {

		// Set up the solution and source arrays
		for(ij=j=0,y=ay;j<n;j++,y+=msu.dy) {
			for(i=0,x=ax;i<m;i++,x+=msu.dx,ij++) {
				z[ij]=0;
				b[ij]=msu.edge(i,ij)?0:1/(x*x+y*y+4);
			}
		}

		// Solve using multigrid V-cycles
		mg.solve_v_cycle(tp);
	}

	// Output the solutions in a format that can be read by Gnuplot using
	// the command "splot 'filename' matrix binary"
	//mg.output_b("b.0");
	//mg.output_z("z.0");
	//mg.output_res("r.0");

	// Delete dynamically allocated memory
	delete [] z;
	delete [] b;
}
