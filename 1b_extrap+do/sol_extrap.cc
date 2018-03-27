#include "sol_extrap.hh"

#include <cstdio>
#include <cstring>

/** Initializes the fourth-order Runge-Kutta solver, allocating memory and
 * setting constants.
 * \param[in] gragg whether to enable the Gragg method.
 * \param[in] dof_ the number of degrees of freedom.
 * \param[in] max_j the maximum size of the extrapolation parameter. */
extrap::extrap(bool gragg_,int dof_,int max_j_) : gragg(gragg_), dof(dof_),
    fcount(0), max_j(max_j_), t(0.), q(new double[dof]), dq(new double[dof*max_j]),
    k1(new double[dof]), aq(gragg?new double[dof]:NULL), n_seq(new int[max_j]) {}

/** The class destructor frees the dynamically allocated memory. */
extrap::~extrap() {
    delete [] n_seq;
    if(gragg) delete [] aq;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

/** Initializes the Romberg sequence for use with the extrapolation routine.*/
void extrap::init_romberg() {
    *n_seq=1;
    for(int j=1;j<max_j;j++) n_seq[j]=n_seq[j-1]*2;
}

/** Initializes the Bulirsch sequence for use with the extrapolation routine.
 */
void extrap::init_bulirsch() {
    int j=1;
    *n_seq=1;
    for(;j<max_j-1;j+=2) {
        n_seq[j]=n_seq[j-1]*2;
        n_seq[j+1]=n_seq[j-1]+n_seq[j];
    }
    if(j<max_j) n_seq[j]=n_seq[j-1]*2;
}

/** Initializes the harmonic sequence for use with the extrapolation routine.
 */
void extrap::init_harmonic() {
    for(int j=0;j<max_j;j++) n_seq[j]=j+1;
}

/** Prints the current state of the solution. */
void extrap::print(double t_,double *in) {
    printf("%g",t_);
    for(int i=0;i<dof;i++) printf(" %g",in[i]);
    puts("");
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step.
 * \param[in] (j,k) the member of the extrapolation family T_{j,k} to use. */
void extrap::solve_fixed(double duration,int steps,bool output,int j,int k) {
    int i,e,l,o;

    // Offset the given j and k to convert to more C++-friendly ranges
    --j;--k;

    // Set up initial condition, compute timestep, and do output if needed
    init();
    double dt=duration/steps;
    if(output) print(t,q);

    // Perform integration steps
    for(l=0;l<steps;l++) {

        // Compute the T_{i,1} terms
        for(i=j-k;i<=j;i++) compute_basic(i,dt);

        // Apply Aitken-Neville algorithm
        for(e=1;e<=k;e++) {
            for(i=j;i>=j-k+e;i--) {
                double *bq=dq+(i-1)*dof,
                       *cq=dq+i*dof,
                       nj=double(n_seq[i]),
                       njk=double(n_seq[i-e]);

                // Square the nj and njk factors if the Gragg method is being
                // used
                if(gragg) {nj*=nj;njk*=njk;}
                double fac=njk/(nj-njk);

                // Use A-N recursion formula
                for(o=0;o<dof;o++) cq[o]=cq[o]+fac*(cq[o]-bq[o]);
            }
        }

        // Copy the required solution into the q array and do output
        // if needed
        memcpy(q,dq+j*dof,dof*sizeof(double));
        t+=dt;
        if(output) print(t,q);
    }
}

/** Computes a basic solution T_{j,1}.
 * \param[in] j the index of timestep subdivision number.
 * \param[in] dt the timestep */
void extrap::compute_basic(int j,double dt) {
    double *cq=dq+j*dof;
    int n=n_seq[j];

    // Set up initial condition for solution
    memcpy(cq,q,dof*sizeof(double));

    // Perform the timesteps
    double t_=t,ldt=dt/n;
    if(gragg) {

        // Perform the first step Gragg update
        ff(t_,cq,k1);
        for(int o=0;o<dof;o++) aq[o]=cq[o]+0.5*ldt*k1[o];
        t_+=ldt;

        // Perform the leapfrogging Gragg updates
        for(int l=1;l<n;l++,t_+=ldt) {

            ff(t_-0.5*ldt,aq,k1);
            for(int o=0;o<dof;o++) cq[o]+=ldt*k1[o];
            ff(t_,cq,k1);
            for(int o=0;o<dof;o++) aq[o]+=ldt*k1[o];
        }

        // Perform final step
        ff(t_-0.5*ldt,aq,k1);
        for(int o=0;o<dof;o++) cq[o]+=ldt*k1[o];
        ff(t_,cq,k1);

        // Assemble smoothed solution S_h(x)
        for(int o=0;o<dof;o++) cq[o]=0.5*(cq[o]+aq[o]+0.5*ldt*k1[o]);
        fcount+=2*n+1;
    } else {

        // Perform the explicit Euler update
        for(int l=0;l<n;l++,t_+=ldt) {
            ff(t_,cq,k1);
            for(int o=0;o<dof;o++) cq[o]+=ldt*k1[o];
        }
        fcount+=n;
    }
}
