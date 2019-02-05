#include <functional>
#include "poisson_fft.hh"
#include "conj_grad.hh"

/** Solve Poisson on [-1,1]x[0,1] by decomposing it into the two glued squares
 *  (1) [-1,0]x[0,1] and (2) [0,1]x[0,1].
 *
 *  The grid looks like:
 *
 *  1 + - - - + - - - +
 *    | . . . | . . . |
 *    | . . . | . . . |
 *    | . . . | . . . |
 *  0 + - - - + - - - +
 *   -1       0       1
 */
class schur : conj_grad
{
    public:
        schur(int n_);
        ~schur();
        void solve(const std::function<double(double,double)>& f);
        virtual void mul_A(double* in, double* out);
        void print_solution();
        /** The number of interior gridpoints per subdomain in one dimension */
        const int n;
        /** The total number of interior gridpoints per subdomain */
        const int nn;
        /** The number of gridpoints on the glue */
        const int ng;
        /** The grid spacing */
        const double h;
        /** The inverse grid spacing squared */
        const double ih2;
    private:
        /** The solver for subdomain 1 */
        poisson_fft grid1;
        /** The solver for subdomain 2 */
        poisson_fft grid2;
        /** The source term on subdomain 1 */
        double* f1;
        /** The source term on subdomain 2 */
        double* f2;
        /** The source term on the glue */
        double* fg;
};
