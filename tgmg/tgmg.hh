#ifndef TGMGPP_HH
#define TGMGPP_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <limits>
#ifdef _OPENMP
#include "omp.h"
#endif

#include "tgmg_config.hh"
#include "tgmg_predict.hh"

// Conversion and output routines for standard types
inline double mod_sq(float c) {return c*c;}
inline float float_output(float c) {return c;}
inline double mod_sq(double c) {return c*c;}
inline float float_output(double c) {return float(c);}
inline double mod_sq(std::complex<double> c) {return std::norm(c);}
inline float float_output(std::complex<double> c) {return float(std::abs(c));}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
inline void tgmg_fatal_error(const char *p,int status) {
    fprintf(stderr,"tgmg: %s\n",p);
    exit(status);
}

/** \brief A general function for setting the accuracy threshold on the squared
 * L2 norm of a multigrid problem, based on a typical scale in the linear
 * system (e.g. the size of a diagonal entry), plus a multiplier (e.g. 1e4)
 * applied to the machine epsilon.
 * \param[in] typical a typical scale.
 * \param[in] mult the multiplier.
 * \return The squared L2 error threshold. */
template<class F>
inline F tgmg_accuracy(F typical,F mult=1e4) {
    F acc=std::numeric_limits<F>::epsilon()*mult*typical;
    return acc*acc;
}

/** \brief Template representing a level of a multigrid hierarchy. */
template<class S,class V,class M>
class tgmg_base {
    public:
        /** The number of grid points in the horizontal direction. */
        const int m;
        /** The number of grid points in the vertical direction. */
        const int n;
        /** The total number of grid points. */
        const int mn;
        /** The horizontal grid points of the child grid (if it
         * exists). */
        const int sm;
        /** The vertical grid points of the child grid (if it exists).
         */
        const int sn;
        /** The periodicity in the x direction. */
        const bool x_prd;
        /** The periodicity in the y direction. */
        const bool y_prd;
        /** The reciprocal of the total number of gridpoints, used in error
         * calculations. */
        const double mn_inv;
        /** A pointer to the source term array. */
        V* b;
        /** A pointer to the solution array. */
        V* z;
        /** A pointer to the source term array for the child grid (if it
         * exists). */
        V* c;
        /** The maximum number of threads used for computations on this level.
         * It will usually be given by the total threads available, but on small
         * grids it may be restricted to ensure each thread has at least two lines
         * to work with. */
        int nthr_max;
        /** The number of threads for smoothing routines. */
        int nthr_smooth;
        /** The number of threads for simple tasks such as clearing the solution. */
        int nthr_simp;
        /** The number of threads for restriction. */
        int nthr_r;
        /** The number of threads for interpolation. */
        int nthr_t;
#ifdef _OPENMP
        /** The number of OpenMP locks. */
        int num_l;
        /** Locks to ensure that Gauss-Seidel is applied consistently. */
        omp_lock_t* ol;
        /** The class destructor destroys the OpenMP locks. */
        ~tgmg_base() {
            if(ol!=NULL) {
                for(int k=0;k<num_l;k++) omp_destroy_lock(ol+k);
                delete [] ol;
            }
        }
        inline void lock(int t) {omp_set_lock(ol+t);}
        inline void unlock(int t) {omp_unset_lock(ol+t);}
#else
        inline void lock(int t) {}
        inline void unlock(int t) {}
#endif
        tgmg_base(S &q_,const int m_,const int n_,const bool x_prd_,const bool y_prd_,V* b_,V* z_)
            : m(m_), n(n_), mn(m*n), sm((m_+(x_prd_?1:2))>>1),
            sn((n_+(y_prd_?1:2))>>1), x_prd(x_prd_), y_prd(y_prd_), mn_inv(1./mn),
            b(b_), z(z_), q(q_) {
#ifdef _OPENMP

            // If OpenMP is available, do a heuristic calculation to guess an
            // appropriate number of threads to use on this level
            nthr_max=omp_get_max_threads();
            if(nthr_max>(n>>1)) nthr_max=n>>1;
            nthr_smooth=thread_trunc(mn>>9);
            nthr_r=thread_trunc(mn>>11);
            nthr_t=thread_trunc(mn>>13);
            nthr_simp=thread_trunc(mn>>15);

            // If the stripped Gauss-Seidel routine is being used, then
            // initialize the OpenMP locks
            if(q.gs_mode<2) ol=NULL;
            else {
                num_l=(y_prd?nthr_max:nthr_max-1);
                ol=new omp_lock_t[num_l];
                for(int k=0;k<num_l;k++) {
                    omp_init_lock(ol+k);
                    omp_set_lock(ol+k);
                }
            }
#else
            // If OpenMP is not available, set the number of threads to 1
            nthr_max=nthr_smooth=nthr_r=nthr_t=nthr_simp=1;
#endif
        }
        void jacobi();
        void zero_jacobi();
        void gauss_seidel();
        void gauss_seidel_reverse();
        void gauss_seidel_from_zero();
        void sor(double omega);
        inline double l2_error() {return mds()*mn_inv;}
        double mds();
        void apply_r();
        void rat();
        void b_to_z();
        void clear_z();
        /** Outputs the solution field to a file.
         * \param[in] filename the name of the file to write to. */
        inline void output_z(const char *filename) {
            output(filename,z,0,1,0,1);
        }
        /** Outputs the source field to a file.
         * \param[in] filename the name of the file to write to. */
        inline void output_b(const char *filename) {
            output(filename,b,0,1,0,1);
        }
        inline void output_res(const char *filename) {
            output_res(filename,0,1,0,1);
        }
        /** Outputs the solution field to a file.
         * \param[in] filename the name of the file to write to. */
        inline void output_z(const char *filename,double ax,double dx,double ay,double dy) {
            output(filename,z,ax,dx,ay,dy);
        }
        /** Outputs the source field to a file.
         * \param[in] filename the name of the file to write to. */
        inline void output_b(const char *filename,double ax,double dx,double ay,double dy) {
            output(filename,b,ax,dx,ay,dy);
        }
        void output_res(const char *filename,double ax,double dx,double ay,double dy);
        void rat_collapse_x();
        void rat_collapse_y();
    protected:
        /** A reference to the multigrid setup class. */
        S &q;
        /** A pointer to the matrix storage on the child grid (if it
         * exists). */
        M* t;
        /** Calculates the residual at a given grid point.
         * \param[in] i the horizontal co-ordinate of the grid point.
         * \param[in] ij the grid point index.
         * \return The residual. */
        inline V res(int i,int ij) {
            return b[ij]-q.mul_a(i,ij)-q.a_cc(i,ij)*z[ij];
        }
    private:
        inline int thread_trunc(int v) {
            return v>=1?(v<=nthr_max?v:nthr_max):1;
        }
        void output(const char *filename,V *ff,double ax,double dx,double ay,double dy);
        void r_double_line_set(V* cp,int ij);
        void r_double_line_add(V* cp,int ij);
        void r_line_set(V* cp,int ij);
        void r_line_add(V* cp,int ij);
        void r_periodic_line_add(V* cp,int ij);
        void rat_line(M* tp,unsigned int wy,int ij);
        void compute_rat(M*& tp,int i,int ij);
        void compute_rat_bdry(M*& tp,unsigned int w,int i,int ij);
};

/** \brief Template representing the lower levels of a multigrid hierarchy. */
template<class V,class M>
class tgmg_level : public tgmg_base<tgmg_level<V,M>,V,M> {
    public:
        using tgmg_base<tgmg_level<V,M>,V,M>::m;
        using tgmg_base<tgmg_level<V,M>,V,M>::n;
        using tgmg_base<tgmg_level<V,M>,V,M>::mn;
        using tgmg_base<tgmg_level<V,M>,V,M>::z;
        using tgmg_base<tgmg_level<V,M>,V,M>::b;
        using tgmg_base<tgmg_level<V,M>,V,M>::x_prd;
        using tgmg_base<tgmg_level<V,M>,V,M>::y_prd;
        using tgmg_base<tgmg_level<V,M>,V,M>::q;
        using tgmg_base<tgmg_level<V,M>,V,M>::t;
        using tgmg_base<tgmg_level<V,M>,V,M>::nthr_max;
        using tgmg_base<tgmg_level<V,M>,V,M>::nthr_t;
        using tgmg_base<tgmg_level<V,M>,V,M>::clear_z;
        using tgmg_base<tgmg_level<V,M>,V,M>::gauss_seidel;
        using tgmg_base<tgmg_level<V,M>,V,M>::gauss_seidel_reverse;
        using tgmg_base<tgmg_level<V,M>,V,M>::gauss_seidel_from_zero;
        /** The calculated matrix entries for the grid. */
        M* const s;
        /** A pointer to the solution array on the parent grid. */
        V* y;
        /** The number of horizontal grid points of the parent grid. */
        const int um;
        /** The number of vertical grid points of the parent grid. */
        const int un;
        /** The mode to use for the Gauss-Seidel smoothing. (0=default) */
        const char gs_mode;
        tgmg_level<V,M>(int m_,int n_,bool x_prd_,bool y_prd_,V* y_,int um_,int un_)
            : tgmg_base<tgmg_level<V,M>,V,M>(*this,m_,n_,x_prd_,y_prd_,new V[m_*n_],new V[m_*n_]),
            s(new M[10*m_*n_]), y(y_), um(um_), un(un_), gs_mode(q.gs_mode) {}
        /** The class destructor clears the dynamically allocated
         * arrays for the solution, source terms, and matrix entries on
         * this level. */
        ~tgmg_level<V,M>() {
            delete [] s;
            delete [] z;
            delete [] b;
        }
        void apply_t();
        inline M a_dl(int i,int ij) {return s[10*ij];}
        inline M a_dc(int i,int ij) {return s[10*ij+1];}
        inline M a_dr(int i,int ij) {return s[10*ij+2];}
        inline M a_cl(int i,int ij) {return s[10*ij+3];}
        inline M a_cc(int i,int ij) {return s[10*ij+4];}
        inline M a_cr(int i,int ij) {return s[10*ij+5];}
        inline M a_ul(int i,int ij) {return s[10*ij+6];}
        inline M a_uc(int i,int ij) {return s[10*ij+7];}
        inline M a_ur(int i,int ij) {return s[10*ij+8];}
        inline V inv_cc(int i,int ij,V v) {return s[10*ij+9]*v;}
        V mul_a(int i,int ij);
        /** Applies Gauss--Seidel sweeps during the first part of the
         * V-cycle, when it is necessary to also set the solution array
         * to zero. If the number of sweeps is non-zero, the rapid
         * zero_jacobi routine is applied in place of the first sweep.
         * \param[in] cyc the number of sweeps to apply. */
        inline void down_gs_iterations(const int cyc) {
            if(cyc==0) clear_z();
            else {
                gauss_seidel_from_zero();
                for(int i=1;i<cyc;i++) gauss_seidel();
            }
        }
    private:
        void t_line(V *yp,V *zp);
};

/** \brief Template representing the top level of a multigrid hierarchy. */
template<class S,class V,class M>
class tgmg : public tgmg_base<S,V,M> {
    public:
        using tgmg_base<S,V,M>::m;
        using tgmg_base<S,V,M>::n;
        using tgmg_base<S,V,M>::mn;
        using tgmg_base<S,V,M>::sm;
        using tgmg_base<S,V,M>::sn;
        using tgmg_base<S,V,M>::rat;
        using tgmg_base<S,V,M>::x_prd;
        using tgmg_base<S,V,M>::y_prd;
        using tgmg_base<S,V,M>::q;
        using tgmg_base<S,V,M>::z;
        using tgmg_base<S,V,M>::c;
        using tgmg_base<S,V,M>::t;
        using tgmg_base<S,V,M>::l2_error;
        using tgmg_base<S,V,M>::jacobi;
        using tgmg_base<S,V,M>::clear_z;
        using tgmg_base<S,V,M>::gauss_seidel;
        using tgmg_base<S,V,M>::gauss_seidel_reverse;
        using tgmg_base<S,V,M>::sor;
        using tgmg_base<S,V,M>::zero_jacobi;
        using tgmg_base<S,V,M>::apply_r;
        using tgmg_base<S,V,M>::nthr_max;
        using tgmg_base<S,V,M>::nthr_smooth;
        using tgmg_base<S,V,M>::nthr_t;
        using tgmg_base<S,V,M>::nthr_r;
        using tgmg_base<S,V,M>::nthr_simp;
        /** The number of child grids in the multigrid hierarchy. */
        int ml;
        /** The verbosity level for status messages. */
        int verbose;
        /** The convergence rate (in digits per iteration) of the previous solve. */
        double conv_rate;
        /** The time for the next thread tuning. */
        double next_tune_time;
        /** An array of pointers to the child grids in the multigrid hierarchy.
         */
        tgmg_level<V,M>* mg[tgmg_max_levels];
        tgmg(S &q_,V* b_,V* z_);
        ~tgmg();
        void print_hierarchy();
        /** Sets up the multigrid hierarchy, computing the matrix entries on
         * all the grids, and performing tuning of the thread numbers on each
         * grid for each operation if needed. */
        inline void setup() {
            calculate_rat();
            try_tuning();
        }
        /** Tries tuning the thread numbers. It will only perform anything if
         * this has not been done since a specified time interval ago. */
        inline void try_tuning() {
#ifdef _OPENMP
            double t0=omp_get_wtime();
            if(t0>next_tune_time) {
                next_tune_time=t0+tgmg_retune_interval;
                tune_threads();
                if(verbose>=1) print_hierarchy();
            }
#endif
        }
        /** Sets up the matrix entries on all grids by recursively conjugating
         * with the restriction and interpolation operators. */
        inline void calculate_rat() {
            rat();for(int l=0;l<ml-1;l++) mg[l]->rat();
        }
        void tune_threads(double res_factor=tgmg_default_resource_factor);
        void thread_speed_test(int l,int type,int &nthr,double lam);
        bool solve(int type,int per_loop,int max_loops,double omega=1,int cyc_down=1,int cyc_up=1,int cyc_bottom=20);
        bool solve(int type,tgmg_predict &tp,double omega=1,int cyc_down=1,int cyc_up=1,int cyc_bottom=20);
        /** Solves the linear system using the Jacobi method. It carries out
         * batches of Jacobi iterations, and checks after each to see if the
         * specified tolerance is reached, after which it terminates.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_jacobi() {
            return solve(0,jacobi_per_loop,max_jacobi_loops);
        }
        /** Solves the linear system using the Gauss--Seidel method. It carries
         * out batches of Gauss--Seidel iterations, and checks after each to
         * see if the specified tolerance is reached, after which it
         * terminates.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_gauss_seidel() {
            return solve(1,gs_per_loop,max_gs_loops);
        }
        /** Solves the linear system using the successive over-relaxation (SOR)
         * method. It carries out batches of SOR iterations, and checks after
         * each to see if the specified tolerance is reached, after which it
         * terminates.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_sor(double omega) {
            return solve(2,sor_per_loop,max_sor_loops,omega);
        }
        /** Solves the linear system using the full multigrid method. It
         * carries out batches of FMG cycles, and checks after each to see if
         * the specified tolerance is reached, after which it terminates.
         * \param[in] (cyc_down,cyc_up,cyc_bottom)
         *            configuration parameters to pass to the v_cycle routine.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_full_multigrid(int cyc_down=1,int cyc_up=1,int cyc_bottom=20) {
            return solve(3,fmg_per_loop,max_fmg_loops,1,cyc_down,cyc_up,cyc_bottom);
        }
        /** Solves the linear system using multigrid V-cycles. It carries out
         * batches of V-cycles, and checks after each to see if the specified
         * tolerance is reached, after which it terminates.
         * \param[in] (cyc_down,cyc_up,cyc_bottom)
         *            configuration parameters to pass to the v_cycle routine.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_v_cycle(int cyc_down=1,int cyc_up=1,int cyc_bottom=20) {
            return solve(4,v_cycle_per_loop,max_v_cycle_loops,1,cyc_down,cyc_up,cyc_bottom);
        }
        /** Solves the linear system using the Jacobi method, using the
         * adaptive approach for choosing iterations.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_jacobi(tgmg_predict &tp) {
            return solve(0,tp);
        }
        /** Solves the linear system using the Gauss--Seidel method, using the
         * adaptive approach for choosing iterations.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_gauss_seidel(tgmg_predict &tp) {
            return solve(1,tp);
        }
        /** Solves the linear system using the successive over-relaxation (SOR)
         * method, using the adaptive approach for choosing iterations.
         * \param[in] omega the over-relaxation parameter.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_sor(tgmg_predict &tp,double omega) {
            return solve(2,tp,omega);
        }
        /** Solves the linear system using the full multigrid method, using
         * the adaptive approach for choosing iterations.
         * \param[in] (cyc_down,cyc_up,cyc_bottom,cyc_top)
         *            configuration parameters to pass to the v_cycle routine.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_full_multigrid(tgmg_predict &tp,int cyc_down=1,int cyc_up=1,int cyc_bottom=20) {
            return solve(3,tp,1,cyc_down,cyc_up,cyc_bottom);
        }
        /** Solves the linear system using multigrid V-cycles, using the
         * adaptive approach for choosing iterations.
         * \param[in] (cyc_down,cyc_up,cyc_bottom,cyc_top)
         *            configuration parameters to pass to the v_cycle routine.
         * \return True if the tolerance was reached before the maximum number
         * of iterations, false otherwise. */
        inline bool solve_v_cycle(tgmg_predict &tp,int cyc_down=1,int cyc_up=1,int cyc_bottom=20) {
            return solve(4,tp,1,cyc_down,cyc_up,cyc_bottom);
        }
        void v_cycle(int cyc_down=1,int cyc_up=1,int cyc_bottom=20);
        void fmg_pre(int cyc_down=1,int cyc_up=1,int cyc_bottom=20);
        void check_symmetric();
        void check_consistent();
    private:
        bool sym_check(int ij,int ij2,M a,M a2,const char* p,const char* p2);
        bool cons_check(int i,int ij,int ij2,M a,const char* p);
        void v_cycle_internal(int level,int cyc_down,int cyc_up,int cyc_bottom);
        void set_threads(int l,int type,int nthr);
        void do_thread_trial(int l,int type);
        inline double iters_and_error(int type,int iters,double omega,int cyc_down,int cyc_up,int cyc_bottom) {
            for(int i=0;i<iters;i++) switch(type) {
                case 0: jacobi();break;
                case 1: gauss_seidel();break;
                case 2: sor(omega);break;
                case 3: fmg_pre(cyc_down,cyc_up,cyc_bottom);
                case 4: v_cycle(cyc_down,cyc_up,cyc_bottom);
            }
            return l2_error();
        }
        inline void iters(int type,int iters,double omega,int cyc_down,int cyc_up,int cyc_bottom) {
            for(int i=0;i<iters;i++) switch(type) {
                case 0: jacobi();break;
                case 1: gauss_seidel();break;
                case 2: sor(omega);break;
                case 3: fmg_pre(cyc_down,cyc_up,cyc_bottom);
                case 4: v_cycle(cyc_down,cyc_up,cyc_bottom);
            }
        }
        inline void iter_message(int i,double acc) {
            printf("Iteration %d, residual %g\n",i,acc);
        }
        inline void bail_message(int i,double iacc,double acc) {
            conv_rate=(log10(iacc)-log10(acc))/i;
            if(verbose==1) printf("Residual=%g not reached %g threshold after %d iterations\n",acc,q.acc,i);
            else if(verbose>=2) printf("%d iters, res %g->%g, %g digits per iter (bailed)\n",i,iacc,acc,conv_rate);
        }
        inline void status_message(int i,double iacc,double acc) {
            conv_rate=(log10(iacc)-log10(acc))/i;
            if(verbose>=2) printf("%d iters, res %g->%g, %g digits per iter\n",i,iacc,acc,conv_rate);
        }
};

#endif
