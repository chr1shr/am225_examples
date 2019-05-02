#include <cstring>

#include "tgmg.hh"

/** The multigrid constructor determines constants from the multigrid setup
 * class, and creates the hierarchy of coarser grids.
 * \param[in] q_ the multigrid setup class to use.
 * \param[in] b_ a pointer to the source array.
 * \param[in] z_ a pointer to the solution array. */
template<class S,class V,class M>
tgmg<S,V,M>::tgmg(S &q_,V* b_,V* z_) : tgmg_base<S,V,M>(q_,q_.m,q_.n,q_.x_prd,q_.y_prd,b_,z_), ml(0),
    verbose(1), conv_rate(std::numeric_limits<double>::max()), next_tune_time(-1.) {
    int am=sm,an=sn,um=m,un=n;
    V* y=z_;

    // Set up the multigrid hierarchy
    while(am*an>tgmg_grid_min) {
        if(ml==tgmg_max_levels) {
            fputs("Maximum levels exceeded\n",stderr);
            exit(TGMGPP_ERROR);
        }
        mg[ml]=new tgmg_level<V,M>(am,an,x_prd,y_prd,y,um,un);
        um=am;am=mg[ml]->sm;
        un=an;an=mg[ml]->sn;
        (ml==0?c:mg[ml-1]->c)=mg[ml]->b;
        (ml==0?t:mg[ml-1]->t)=mg[ml]->s;
        y=mg[ml]->z;
        ml++;
    }
}

/** The multigrid destructor frees the memory used for the grid hierarchy. */
template<class S,class V,class M>
tgmg<S,V,M>::~tgmg() {
    while(ml>0) delete mg[--ml];
}

/** Prints the grid hierarchy. */
template<class S,class V,class M>
void tgmg<S,V,M>::print_hierarchy() {
    printf("# Top grid level: (%d,%d) [%s,%s] {%d,%d,-,%d}\n",
           m,n,m&1?"odd":"even",n&1?"odd":"even",nthr_smooth,nthr_r,nthr_simp);
    for(int l=0;l<ml;l++) {
        int gm=mg[l]->m,gn=mg[l]->n;
        printf("# Grid level %2d : (%d,%d) [%s,%s] {%d,",
               l,gm,gn,gm&1?"odd":"even",gn&1?"odd":"even",mg[l]->nthr_smooth);
        if(l==ml-1) printf("-,%d,%d}\n",mg[l]->nthr_t,mg[l]->nthr_simp);
        else printf("%d,%d,%d}\n",mg[l]->nthr_r,mg[l]->nthr_t,mg[l]->nthr_simp);
    }
}

/** Calculates the sum of squares of residuals.
 * \return The sum. */
template<class S,class V,class M>
double tgmg_base<S,V,M>::mds() {
    double c=0;
#pragma omp parallel for num_threads(nthr_r) reduction(+:c)
    for(int j=0;j<mn;j+=m) {
        double cl=0;
        for(int i=0,ij=j;i<m;i++,ij++) cl+=mod_sq(b[ij]-q.mul_a(i,ij)-q.a_cc(i,ij)*z[ij]);
        c+=cl;
    }
    return c;
}

/** Carries out a Jacobi iteration, by creating a temporary array to store the
 * improved solution, and then copying it over the main solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::jacobi() {
    V *zn=new V[mn];
#pragma omp parallel num_threads(nthr_smooth)
    {
#pragma omp for
        for(int j=0;j<mn;j+=m) {
            for(int i=0,ij=j;i<m;i++,ij++) zn[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
        }
#pragma omp for
        for(int j=0;j<mn;j+=m) {
            for(V *znp=zn+j,*zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=*(znp++);
        }
    }
    delete [] zn;
}

/** Carries out a Gauss-Seidel sweep, assuming that the solution was zero beforehand. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::gauss_seidel_from_zero() {
    if(q.gs_mode<2) {
#pragma omp parallel num_threads(nthr_smooth)
       {
#pragma omp for
            for(int j=m;j<mn;j+=2*m) {
                z[j]=q.inv_cc(0,j,b[j]);
                int ij=j+1;
                for(int i=1;i<m-1;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.a_cl(i,ij)*z[ij-1]);
                z[ij]=q.inv_cc(m-1,ij,(x_prd?b[ij]-q.a_cr(m-1,ij)*z[ij-m+1]:b[ij])-q.a_cl(m-1,ij)*z[ij-1]);
            }
            if(y_prd&&(n&1)) {
                for(V* zp=z;zp<z+m;zp++) *zp=V(0.);
                for(V* zp=z+mn-m;zp<z+mn;zp++) *zp=V(0.);
            }
#pragma omp for
            for(int j=0;j<mn;j+=2*m) {
                int ij=j;
                if(x_prd) z[ij+m-1]=V(0.);
                for(int i=0;i<m-1;i++,ij++) {z[ij+1]=0.;z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));}
                z[ij]=q.inv_cc(m-1,ij,b[ij]-q.mul_a(m-1,ij));
            }
       }
    } else {
        clear_z();
        gauss_seidel();
    }
}

/** Carries out a Gauss-Seidel sweep. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::gauss_seidel() {
#pragma omp parallel num_threads(nthr_smooth)
    {
        if(q.gs_mode<2) {
#pragma omp for
            for(int j=m;j<mn;j+=2*m) {
                for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
            }
#pragma omp for
            for(int j=0;j<mn;j+=2*m) {
                for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
            }
        } else {

            // Stripped Gauss-Seidel, using OpenMP locks to ensure that the
            // strips don't conflict with each other
            int t=omp_get_thread_num(),js=t*n/nthr_smooth*m,je=(t+1)*n/nthr_smooth*m;
            for(int j=js;j<je;j+=m) {
                if(j==je-m&&(y_prd||t<nthr_smooth-1)) lock(t);
                for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
                if(j==js&&(y_prd||t!=0)) unlock(t==0?nthr_smooth-1:t-1);
            }
        }
    }
}

/** Carries out a Gauss-Seidel sweep in reverse on the grid */
template<class S,class V,class M>
void tgmg_base<S,V,M>::gauss_seidel_reverse() {
#pragma omp parallel num_threads(nthr_smooth)
    {
        if(q.gs_mode<2) {
#pragma omp for
            for(int j=2*m-1;j<mn;j+=2*m) {
                for(int i=m-1,ij=j;i>=0;i--,ij--) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
            }
#pragma omp for
            for(int j=m-1;j<mn;j+=2*m) {
                for(int i=m-1,ij=j;i>=0;i--,ij--) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
            }
        } else {
            int t=omp_get_thread_num(),js=t*n/nthr_smooth*m,je=(t+1)*n/nthr_smooth*m;
            for(int j=je-1;j>=js;j-=m) {
                if(j==js+m-1&&(y_prd||t!=0)) lock(num_l+(t==0?nthr_smooth-1:t-1));
                for(int i=m-1,ij=j;i>=0;i--,ij--) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
                if(j==je-1&&(y_prd||t<nthr_smooth-1)) unlock(num_l+t);
            }
        }
    }
}

/** Carries out a successive over-relation (SOR) sweep.
 * \param[in] omega the factor by which to relax by, usually between 0 and 2.
 */
template<class S,class V,class M>
void tgmg_base<S,V,M>::sor(double omega) {
#pragma omp parallel num_threads(nthr_smooth)
    {
#pragma omp for
        for(int j=0;j<mn;j+=2*m) {
            for(int i=0,ij=j;i<m;i++,ij++) z[ij]+=omega*(q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij))-z[ij]);
        }
#pragma omp for
        for(int j=m;j<mn;j+=2*m) {
            for(int i=0,ij=j;i<m;i++,ij++) z[ij]+=omega*(q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij))-z[ij]);
        }
    }
}

/** Solves the linear system using one of the smoothing techniques. It performs
 * batches of smoothing cycles and checks after each to see if the specified
 * tolerance is reached, after which it terminates.
 * \param[in] type the type of smoothing to perform.
 * \param[in] per_loop the number of smoothing cycles to perform in each batch.
 * \param[in] omega the over-relaxation parameter, only used with SOR smoothing.
 * \param[in] (cyc_down,cyc_up,cyc_bottom) V-cycle configuration parameters.
 * \return True if the tolerance was reached before the maximum number of
 * iterations, false otherwise. */
template<class S,class V,class M>
bool tgmg<S,V,M>::solve(int type,int per_loop,int max_loops,double omega,int cyc_down,int cyc_up,int cyc_bottom) {
    double iacc=1,acc=0;int n=0;
    if(verbose>=2) iacc=l2_error();
    do {
        if(++n==max_loops) {
            bail_message(n*per_loop,iacc,acc);
            return false;
        }
        acc=iters_and_error(type,per_loop,omega,cyc_down,cyc_up,cyc_bottom);
        if(verbose==3) iter_message(n*per_loop,acc);
    } while(acc>q.acc);
    status_message(n*per_loop,iacc,acc);
    return !std::isnan(acc);
}

/** Solves the linear system using one of the smoothing techniques and the
 * adaptive approach for choosing iterations.
 * \param[in] type the type of smoothing to perform.
 * \param[in] tp a class for predicting the number of smoothing steps needed.
 * \param[in] omega the over-relaxation parameter, only used with SOR smoothing.
 * \param[in] (cyc_down,cyc_up,cyc_bottom) V-cycle configuration parameters.
 * \return True if the tolerance was reached before the maximum number of
 * iterations, false otherwise. */
template<class S,class V,class M>
bool tgmg<S,V,M>::solve(int type,tgmg_predict &tp,double omega,int cyc_down,int cyc_up,int cyc_bottom) {
    int n=tp.lim/tp.mult;

    // Perform the predicted number of smoothing iterations
    double iacc=verbose>=2?l2_error():1,
           acc=iters_and_error(type,n,omega,cyc_down,cyc_up,cyc_bottom);
    if(verbose==3) iter_message(n,acc);

    // Test the L2 error against the given threshold. If it's lower, then try
    // to decrease the smoothing iterations for the solve.
    if(acc<q.acc) {
        if(tp.lim>0) tp.lim-=1+tp.lim/tp.decay;
    } else {

        // If not, then try more smoothing operations, testing the
        // accuracy after every triangular number of iterations
        int k=0;
        do {
            tp.lim+=++k*tp.mult;
            if(tp.lim>tp.max_thresh) {
                bail_message(n,iacc,acc);
                return false;
            }
            n+=k;
            acc=iters_and_error(type,k,omega,cyc_down,cyc_up,cyc_bottom);
            if(verbose==3) iter_message(n,acc);
        } while(acc>=q.acc);
    }

    // Record the number of iterations and perform any extra smoothing
    // steps
    tp.add_iters(n);
    status_message(n,iacc,acc);
    iters(type,tp.extra_iters,omega,cyc_down,cyc_up,cyc_bottom);
    return !std::isnan(acc);
}

/** Carries out a V-cycle.
 * \param[in] cyc_down the number of Gauss--Seidel sweeps to apply on the way
 *                     down the grid hierarchy.
 * \param[in] cyc_up the number of Gauss--Seidel sweeps to apply on the way up
 *                   the grid hierarchy, not including the top level.
 * \param[in] cyc_bottom the number of Gauss--Seidel sweeps to apply on the
 *                       bottom level of the grid hierarchy. */
template<class S,class V,class M>
void tgmg<S,V,M>::v_cycle(int cyc_down,int cyc_up,int cyc_bottom) {
    for(int j=0;j<cyc_down;j++) gauss_seidel();

    // Only do the bulk of the V-cycle if levels are actually defined
    if(ml>0) {

        // Propagate the solution down the hierarchy, smoothing at each step
        apply_r();
        mg[0]->down_gs_iterations(cyc_down);
        v_cycle_internal(0,cyc_down,cyc_up,cyc_bottom);
        mg[0]->apply_t();
    }

    // Apply smoothing sweeps on the top level
    for(int j=0;j<cyc_up;j++) gauss_seidel();
}

/** Performs a partial V-cycle from a given level in the hierarchy.
 * \param[in] level the level to start from.
 * \param[in] cyc_down the number of Gauss--Seidel sweeps to apply on the way
 *                     down the grid hierarchy.
 * \param[in] cyc_up the number of Gauss--Seidel sweeps to apply on the way up
 *                   the grid hierarchy, not including the top level.
 * \param[in] cyc_bottom the number of Gauss--Seidel sweeps to apply on the
 *                       bottom level of the grid hierarchy. */
template<class S,class V,class M>
void tgmg<S,V,M>::v_cycle_internal(int level,int cyc_down,int cyc_up,int cyc_bottom) {
    int i,j;

    // Propagate the solution down the hierarchy, smoothing at each step
    for(i=level;i<ml-1;i++) {
        if(i>level) mg[i]->down_gs_iterations(cyc_down);
        mg[i]->apply_r();
    }

    // Carry out smoothing sweeps on the bottom level
    mg[ml-1]->down_gs_iterations(cyc_bottom);

    // Propagate the solution up the hierarchy, smoothing at each step
    for(i=ml-1;i>level;i--) {
        mg[i]->apply_t();
        for(j=0;j<cyc_up;j++) mg[i-1]->gauss_seidel();
    }
}

/** Carries out a full multigrid iteration, apart from the final V-cycle.
 * \param[in] cyc_down the number of Gauss--Seidel sweeps to apply on the way
 *                     down the grid hierarchy.
 * \param[in] cyc_up the number of Gauss--Seidel sweeps to apply on the way up
 *                   the grid hierarchy, not including the top level.
 * \param[in] cyc_bottom the number of Gauss--Seidel sweeps to apply on the
 *                       bottom level of the grid hierarchy. */
template<class S,class V,class M>
void tgmg<S,V,M>::fmg_pre(int cyc_down,int cyc_up,int cyc_bottom) {
    int i,j;
    for(j=0;j<cyc_down;j++) gauss_seidel();

    // Only do the bulk of the V-cycle if levels are actually defined
    if(ml>0) {

        // Propagate the solution down the hierarchy, smoothing at each step
        apply_r();
        for(i=0;i<ml-1;i++) {
            mg[i]->down_gs_iterations(cyc_down);
            mg[i]->apply_r();
        }

        // Carry out smoothing sweeps on the bottom level
        mg[ml-1]->down_gs_iterations(cyc_bottom);

        // Propagate the solution up the hierarchy, applying a partial V-cycle
        // at each step
        for(i=ml-1;i>0;i--) {
            mg[i]->apply_t();
            for(j=0;j<cyc_up+cyc_down;j++) mg[i-1]->gauss_seidel();
            v_cycle_internal(i-1,cyc_down,cyc_up,cyc_bottom);
        }

        // Interpolate to the top level
        mg[0]->apply_t();
    }

    // Apply a final complete V-cycle
    for(j=0;j<cyc_up;j++) gauss_seidel();
}

/** Tunes the number of threads used in the multigrid hierarchy by performing timing tests.
 * \param[in] res_factor the factor by which to weight total consumed system
 *                       resources along with wall clock time. */
template<class S,class V,class M>
void tgmg<S,V,M>::tune_threads(double res_factor) {
    V* temp_z=new V[mn];
    int nthr[4];
    for(int type=0;type<4;type++) nthr[type]=1;

    for(int l=ml-1;l>=-1;l--) {

        // Move the solution to temporary space
        V* cur_z=l==-1?z:mg[l]->z;
        size_t size_z=(l==-1?mn:mg[l]->mn)*sizeof(V);
        memcpy(temp_z,cur_z,size_z);

        // Tune all four operations that use this solution array
        for(int type=0;type<(l==ml-1?2:4);type++)
            thread_speed_test(l,type,nthr[type],res_factor);

        // Put back the solution
        memcpy(cur_z,temp_z,size_z);
    }
    delete [] temp_z;
}

/** Selects the number of threads to use for particular operation based on speed tests.
 * \param[in] l the level to consider (with -1 being this level).
 * \param[in] type the operation to consider (0: smoothing, 1: simple, 2:
 *                 interpolation, 3:restriction).
 * \param[in,out] nthr the maximum number of threads to test; on output, this
 *                     value is set to the chosen optimal value.
 * \param[in] res_factor the factor by which to weight total consumed system
 *                       resources along with wall clock time. */
template<class S,class V,class M>
void tgmg<S,V,M>::thread_speed_test(int l,int type,int &nthr,double res_factor) {
#ifdef _OPENMP
    int max_thr=type==2?mg[l+1]->nthr_max:(l==-1?nthr_max:mg[l]->nthr_max);
    if(nthr==1&&max_thr==nthr_max) {
        const double large=std::numeric_limits<double>::max();
        double min_f=large,f;

        // From the starting point, test decreasing numbers of threads
        for(int o=0;o<2;o++) {
            int k=o==0?nthr_max:1,s;
            if(o==0) do_thread_trial(l,type);
            set_threads(l,type,k);
            double t0=omp_get_wtime(),tbase=t0,del_t,min_t=large,t1;
            do {
                s=0;

                // Run a group timing tests, up to the specified duration
                do {
                    do_thread_trial(l,type);
                    s++;t1=omp_get_wtime();del_t=t1-t0;
                } while(del_t<tgmg_time_group);

                // Store the minimum time for the group
                t0=t1;
                del_t/=s;
                if(del_t<min_t) min_t=del_t;
            } while(t0-tbase<tgmg_time_max);

            // Compute measure of efficiency, using a combination of the
            // measured time and the total consumed system resources
            f=min_t*(1+res_factor*k);
            if(f<min_f) {min_f=f;nthr=k;}

            // If the efficiency measure is way off the best seen, then bail
            // out and skip the remaining tests
            if(f>tgmg_f_bail_factor*min_f) break;
        }
    }
    set_threads(l,type,nthr);
#endif
}

/** Performs a multigrid operation that operates on the solution vector on a
 * given level, for the purpose of thread tuning.
 * \param[in] l the level to consider (with -1 being this level).
 * \param[in] type the operation to consider (0: smoothing, 1: simple, 2:
 *                 interpolation, 3:restriction). */
template<class S,class V,class M>
void tgmg<S,V,M>::do_thread_trial(int l,int type) {
    if(l==-1) {
        switch(type) {
            case 0: gauss_seidel();break;
            case 1: zero_jacobi();break;
            case 2: (*mg)->apply_t();break;
            case 3: apply_r();
        }
    } else {
        switch(type) {
            case 0: mg[l]->gauss_seidel();break;
            case 1: mg[l]->zero_jacobi();break;
            case 2: mg[l+1]->apply_t();break;
            case 3: mg[l]->apply_r();
        }
    }
}

/** Sets the threads to use when operating on the solution vector on a
 * given level, for the purpose of thread tuning.
 * \param[in] l the level to consider (with -1 being this level).
 * \param[in] type the operation to consider (0: smoothing, 1: simple, 2:
 *                 interpolation, 3:restriction).
 * \param[in] nthr the number of threads to set. */
template<class S,class V,class M>
void tgmg<S,V,M>::set_threads(int l,int type,int nthr) {
    if(l==-1) {
        switch(type) {
            case 0: nthr_smooth=nthr;break;
            case 1: nthr_simp=nthr;break;
            case 2: (*mg)->nthr_t=nthr;break;
            case 3: nthr_r=nthr;
        }
    } else {
        switch(type) {
            case 0: mg[l]->nthr_smooth=nthr;break;
            case 1: mg[l]->nthr_simp=nthr;break;
            case 2: mg[l+1]->nthr_t=nthr;break;
            case 3: mg[l]->nthr_r=nthr;
        }
    }
}

/** Calculates the matrix product f$A'z\f$ at a grid point where \f$A'\f$ is a
 * matrix of all off-diagonal entries of \f$A\f$.
 * \param[in] i the horizontal co-ordinate of the grid point to consider.
 * \param[in] ij the grid point index.
 * \return The matrix product at the grid point. */
template<class V,class M>
V tgmg_level<V,M>::mul_a(int i,int ij) {
    M *e=s+10*ij;V *f=z+ij;
    return ij<m?(i==0?e[5]*f[1]+e[7]*f[m]+e[8]*f[m+1]+(x_prd?e[3]*f[m-1]+e[6]*f[2*m-1]
                   +(y_prd?*e*f[mn-1]:V(0.0)):(y_prd?e[1]*f[mn-m]+e[2]*f[mn-m+1]:V(0.0))):
             i==m-1?e[3]*f[-1]+e[6]*f[m-1]+e[7]*f[m]+(x_prd?e[5]*f[1-m]+e[8]*f[1]
                   +(y_prd?e[2]*z[mn-m]:V(0.0)):(y_prd?*e*f[mn-m-1]+e[1]*f[mn-m]:V(0.0))):
             e[3]*f[-1]+e[5]*f[1]+e[6]*f[m-1]+e[7]*f[m]+e[8]*f[m+1]
                   +(y_prd?*e*f[mn-m-1]+e[1]*f[mn-m]+e[2]*f[mn-m+1]:V(0.0))):
           ij>=mn-m?(i==0?e[1]*f[-m]+e[2]*f[1-m]+e[5]*f[1]+(x_prd?*e*f[-1]+e[3]*f[m-1]
                   +(y_prd?e[6]*z[m-1]:V(0.0)):(y_prd?e[7]*f[m-mn]+e[8]*f[m+1-mn]:V(0.0))):
             i==m-1?*e*f[-m-1]+e[1]*f[-m]+e[3]*f[-1]+(x_prd?e[2]*f[1-2*m]+e[5]*f[1-m]
                   +(y_prd?e[8]*(*z):V(0.0)):(y_prd?e[6]*f[m-1-mn]+e[7]*f[m-mn]:V(0.0))):
             *e*f[-m-1]+e[1]*f[-m]+e[2]*f[1-m]+e[3]*f[-1]+e[5]*f[1]
                   +(y_prd?e[6]*f[m-1-mn]+e[7]*f[m-mn]+e[8]*f[m+1-mn]:V(0.0))):
           (i==0?e[1]*f[-m]+e[2]*f[1-m]+e[5]*f[1]+e[7]*f[m]+e[8]*f[m+1]
                   +(x_prd?*e*f[-1]+e[3]*f[m-1]+e[6]*f[2*m-1]:V(0.0)):
        i==m-1?*e*f[-m-1]+e[1]*f[-m]+e[3]*f[-1]+e[6]*f[m-1]+e[7]*f[m]
                   +(x_prd?e[2]*f[1-2*m]+e[5]*f[1-m]+e[8]*f[1]:V(0.0)):
        *e*f[-m-1]+e[1]*f[-m]+e[2]*f[1-m]+e[3]*f[-1]+e[5]*f[1]
                   +e[6]*f[m-1]+e[7]*f[m]+e[8]*f[m+1]);
}

/** Applies the restriction operation, computing the residual at the current
 * level and projecting it the level below. */
template<class V,class M>
void tgmg_level<V,M>::apply_t() {

    // Deal with the bulk of the grid
#pragma omp parallel for num_threads(nthr_t)
    for(int j=0;j<n-(y_prd||(un&1)?1:2);j++) {
        V *yp=y+((j*um)<<1),*zp=z+j*m,*ze=zp+(m-(x_prd||(um&1)?1:2)),zi,zi2;
        while(zp<ze) {
            *yp+=*zp;zi=*zp+zp[m];
            yp[um]+=0.5*zi;yp++;
            zi2=*(zp++);zi+=*zp+zp[m];
            zi2+=*zp;
            *yp+=0.5*zi2;yp[um]+=0.25*zi;yp++;
        }
        *yp+=*zp;yp[um]+=0.5*(*zp+zp[m]);
        if(!(um&1)) {
            yp++;
            if(x_prd) {
                zi=*zp+zp[1-m];
                *yp+=0.5*zi;yp[um]+=0.25*(zi+zp[m]+zp[1]);
            } else {
                zp++;
                *yp+=*zp;yp[um]+=0.5*(*zp+zp[m]);
            }
        }
    }

    // Deal with the final line of the grid
    if(un&1) t_line(y+um*(un-1),z+m*(n-1));
    else {
        V *yp=y+um*(un-1),*zp=z+m*(n-1);
        if(y_prd) {
            t_line(yp-um,zp);
            V *ze=zp+(m-(x_prd||(um&1)?1:2)),zi;
            while(zp<ze) {
                zi=*zp+zp[m-mn];
                *(yp++)+=0.5*zi;zp++;
                *yp+=0.25*(zi+*zp+zp[m-mn]);yp++;
            }
            *yp+=0.5*(*zp+zp[m-mn]);
            if(!(um&1)) {
                yp++;
                if(x_prd) {
                    zi=*zp+zp[1-m];
                    *yp+=0.25*(*zp+zp[1-m]+zp[m-mn]+*z);
                } else {
                    zp++;
                    *yp+=0.5*(*zp+zp[m-mn]);
                }
            }
        } else {
            t_line(yp-um,zp-m);
            t_line(yp,zp);
        }
    }
}

/** Interpolates a single line from the current grid to the parent grid.
 * \param[in] yp a pointer to the start of the line of the parent grid.
 * \param[in] zp a pointer to the start of the line of the current grid. */
template<class V,class M>
void tgmg_level<V,M>::t_line(V *yp,V* zp) {
    V *ze=zp+(m-(x_prd||(um&1)?1:2)),zi;
    while(zp<ze) {
        *(yp++)+=*zp;
        zi=*(zp++);
        *(yp++)+=0.5*(zi+*zp);
    }
    *yp+=*zp;
    if(!(um&1)) yp[1]+=x_prd?0.5*(*zp+zp[1-m]):zp[1];
}

/** Calculates the residual at the current level and restricts it to the child
 * grid. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::apply_r() {

    // Handle the first bulk pass where the residuals are set into the child
    // array
#pragma omp parallel for num_threads(nthr_r)
    for(int j=0;j<n-2;j+=4) r_double_line_set(c+sm*(j>>1),j*m);

    // Handle the boundary lines, taking into account periodicity if necessary
    V *cp=c+((n-1)>>1)*sm;
    int ij=m*((n-1)&~1);
    if((y_prd||(n&1)?sn:sn-1)&1) r_line_set(cp,ij);else r_line_add(cp,ij);
    if((n&1)==0) {
        ij+=m;
        if(y_prd) r_periodic_line_add(cp,ij);
        else r_line_set(cp+sm,ij);
    }

    // Handle the second bulk pass where the residuals are added to the
    // existing values in the child array
#pragma omp parallel for num_threads(nthr_r)
    for(int j=2;j<n-2;j+=4) r_double_line_add(c+sm*(j>>1),j*m);
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_double_line_set(V* cp,int ij) {
    V ci,ci2;
    ci=0.5*res(0,ij+m);
    cp[sm]=ci;*cp=res(0,ij)+ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
        cp[sm]+=ci2;*(cp++)+=ci+ci2;
        i++;ij++;
        ci2+=0.5*res(i,ij+m);cp[sm]=ci2;
        *cp=res(i,ij)+ci+ci2;
        i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {
            ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
            cp[sm]+=ci2;*cp+=ci+ci2;
            cp[1]+=ci2;cp[1-sm]+=ci+ci2;
        } else {
            cp++;ci=0.5*res(i,ij+m);
            *cp=res(i,ij)+ci;cp[sm]=ci;
        }
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_double_line_add(V* cp,int ij) {
    V ci,ci2;
    ci=0.5*res(0,ij+m);
    cp[sm]+=ci;*cp+=res(0,ij)+ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
        cp[sm]+=ci2;*(cp++)+=ci+ci2;
        i++;ij++;
        ci2+=0.5*res(i,ij+m);cp[sm]+=ci2;
        *cp+=res(i,ij)+ci+ci2;
        i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {
            ci=0.5*res(i,ij);ci2=0.25*res(i,ij+m);
            cp[sm]+=ci2;*cp+=ci+ci2;
            cp[1]+=ci2;cp[1-sm]+=ci+ci2;
        } else {
            cp++;ci=0.5*res(i,ij+m);
            *cp+=res(i,ij)+ci;cp[sm]+=ci;
        }
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_line_set(V* cp,int ij) {
    V ci=res(0,ij);
    *cp=ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);*cp+=ci;cp++;i++;ij++;
        ci+=res(i,ij);*cp=ci;i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {ci=0.5*res(i,ij);*cp+=ci;cp[1-sm]+=ci;}
        else {cp++;ci=res(i,ij);*cp=ci;}
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_line_add(V* cp,int ij) {
    V ci=res(0,ij);
    *cp+=ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.5*res(i,ij);*cp+=ci;cp++;i++;ij++;
        ci+=res(i,ij);*cp+=ci;i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {ci=0.5*res(i,ij);*cp+=ci;cp[1-sm]+=ci;}
        else {cp++;ci=res(i,ij);*cp+=ci;}
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_periodic_line_add(V* cp,int ij) {
    const int sd=sm*(1-sn);
    V ci=0.5*res(0,ij);
    *cp+=ci;cp[sd]+=ci;
    int i=1;ij++;
    while(i<m-1) {
        ci=0.25*res(i,ij);*cp+=ci;cp[sd]+=ci;cp++;i++;ij++;
        ci+=0.5*res(i,ij);*cp+=ci;cp[sd]+=ci;i++;ij++;
    }
    if(i==m-1) {
        if(x_prd) {ci=0.25*res(i,ij);*cp+=ci;cp[sd]+=ci;cp[1-sm]+=ci;*c+=ci;}
        else {cp++;ci=0.5*res(i,ij);*cp+=ci;cp[sd]+=ci;}
    }
}

/** Calculates the matrix entries on the child grid by conjugating the matrix
 * entries on this grid by the restriction and interpolation operators. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::rat() {
    if(y_prd) {
        if((n&1)==0) {
            rat_line(t,176,0); //IR
            if(n>2) {
#pragma omp parallel for num_threads(nthr_smooth)
                for(int j=2;j<n;j+=2) rat_line(t+5*sm*j,160,j*m);
            } else rat_collapse_y();
        } else if(n==1) rat_line(t,0,0); //HH
        else {
            rat_line(t,144,0); //DR
#pragma omp parallel for num_threads(nthr_smooth)
            for(int j=2;j<n-2;j+=2) rat_line(t+5*sm*j,160,j*m); //RR
            rat_line(t+10*sm*(sn-1),96,(n-1)*m); //RD
        }
    } else {
        if(n==2) {
            rat_line(t,64,0); //HD
            rat_line(t+10*sm,16,1); //DH
        } else {
            rat_line(t,128,0); //HR
#pragma omp parallel for num_threads(nthr_smooth)
            for(int j=2;j<n-2;j+=2) rat_line(t+5*sm*j,160,j*m); //RR
            if(n&1) rat_line(t+10*sm*(sn-1),32,(n-1)*m); //RH
            else {
                rat_line(t+10*sm*(sn-2),96,(n-2)*m); //RD
                rat_line(t+10*sm*(sn-1),16,(n-1)*m); //DH
            }
        }
    }
    if(x_prd&&m==2) rat_collapse_x();
}

/** Calculates the elements of the matrix for a horizontal line in the child
 * grid, by conjugating with the restriction and interpolation operators.
 * \param[in] tp the memory location to start storing the matrix elements at.
 * \param[in] wy a mask giving the boundary information in the y direction.
 * \param[in] ij the grid point index of the first point in the line to
 *               consider. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::rat_line(M* tp,unsigned int wy,int ij) {
    int i;
    if(x_prd) {
        if((m&1)==0) {
            compute_rat_bdry(tp,wy|11,0,ij);ij+=2; //IR
            if(wy==160) {for(i=2;i<m;i+=2,ij+=2) compute_rat(tp,i,ij);}
            else for(i=2;i<m;i+=2,ij+=2) compute_rat_bdry(tp,wy|10,i,ij); //RR
        } else if(m==1) compute_rat_bdry(tp,wy,0,ij); //HH
        else {
            compute_rat_bdry(tp,wy|9,0,ij);ij+=2; //DR
            if(wy==160) {for(i=2;i<m-2;i+=2,ij+=2) compute_rat(tp,i,ij);}
            else for(i=2;i<m-2;i+=2,ij+=2) compute_rat_bdry(tp,wy|10,i,ij); //RR
            compute_rat_bdry(tp,wy|6,i,ij); //RD
        }
    } else {
        if(m==2) {
            compute_rat_bdry(tp,wy|4,0,ij); //HD
            compute_rat_bdry(tp,wy|1,1,ij+1); //DH
        } else {
            compute_rat_bdry(tp,wy|8,0,ij);ij+=2; //HR
            if(wy==160) {for(i=2;i<m-2;i+=2,ij+=2) compute_rat(tp,i,ij);}
            else for(i=2;i<m-2;i+=2,ij+=2) compute_rat_bdry(tp,wy|10,i,ij); //RR
            if(m&1) compute_rat_bdry(tp,wy|2,i,ij); //RH
            else {
                compute_rat_bdry(tp,wy|6,i,ij); //RD
                compute_rat_bdry(tp,wy|1,i+1,ij+1); //DH
            }
        }
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::rat_collapse_x() {
    printf("X %d %d %d %d\n",m,n,sm,sn);
    for(M *tp=t;tp<t+10*sn;tp+=10) {
        tp[1]+=*tp+tp[2];*tp=M(0.0);tp[2]=M(0.0);
        tp[4]+=tp[3]+tp[5];tp[9]=1./tp[4];tp[3]=M(0.0);tp[5]=M(0.0);
        tp[7]+=tp[6]+tp[8];tp[6]=M(0.0);tp[8]=M(0.0);
    }
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::rat_collapse_y() {
    printf("Y %d %d %d %d\n",m,n,sm,sn);
    for(M *tp=t;tp<t+10*sm;tp+=10) {
        tp[3]+=*tp+tp[6];*tp=M(0.0);tp[6]=M(0.0);
        tp[4]+=tp[1]+tp[7];tp[9]=1./tp[4];tp[1]=M(0.0);tp[7]=M(0.0);
        tp[5]+=tp[2]+tp[8];tp[2]=M(0.0);tp[8]=M(0.0);
    }
}

/** Calculates the elements of the matrix at one grid point on the next level
 * down the hierarchy by conjugating with the restriction and interpolation
 * operators. This routine is specifically for a grid point in the bulk where
 * boundaries do not need to be accounted for.
 * \param[in,out] tp a reference to memory in which to store the matrix
 *                   elements.
 * \param[in] i the horizontal co-ordinate of the grid point.
 * \param[in] ij the grid point index. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::compute_rat(M*& tp,int i,int ij) {
    int ci,cij;
    M dl,dc,dr,cl,cc,cr,ul,uc,ur;

    // Bottom left contribution
    cij=ij-m-1;ci=i-1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
    cr+=cc;cc+=cl;ur+=uc;uc+=ul;
    *tp=cc;tp[1]=cr;tp[3]=uc;tp[4]=ur;

    // Bottom center contribution
    cij=ij-m;
    dl=0.25*q.a_dl(i,cij);dc=0.5*q.a_dc(i,cij);dr=0.25*q.a_dr(i,cij);
    cl=0.125*q.a_cl(i,cij);cc=0.25*q.a_cc(i,cij);cr=0.125*q.a_cr(i,cij);
    ul=0.25*q.a_ul(i,cij);uc=0.5*q.a_uc(i,cij);ur=0.25*q.a_ur(i,cij);
    ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
    cc+=cl;uc+=ul;
    cc+=cr;uc+=ur;
    *tp+=cl;tp[1]+=cc;tp[2]=cr;tp[3]+=ul;tp[4]+=uc;tp[5]=ur;

    // Bottom right contribution
    cij=ij-m+1;ci=i+1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
    cl+=cc;cc+=cr;ul+=uc;uc+=ur;
    tp[1]+=cl;tp[2]+=cc;tp[4]+=ul;tp[5]+=uc;

    // Middle left contribution
    cij=ij-1;ci=i-1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.5*q.a_cl(ci,cij);cc=0.25*q.a_cc(ci,cij);cr=0.5*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    cl+=dl;cc+=dc;cr+=dr;
    cl+=ul;cc+=uc;cr+=ur;
    dr+=dc;dc+=dl;cr+=cc;cc+=cl;ur+=uc;uc+=ul;
    *tp+=dc;tp[1]+=dr;tp[3]+=cc;tp[4]+=cr;tp[6]=uc;tp[7]=ur;

    // Middle center contribution
    dl=0.25*q.a_dl(i,ij);dc=0.5*q.a_dc(i,ij);dr=0.25*q.a_dr(i,ij);
    cl=0.5*q.a_cl(i,ij);cc=q.a_cc(i,ij);cr=0.5*q.a_cr(i,ij);
    ul=0.25*q.a_ul(i,ij);uc=0.5*q.a_uc(i,ij);ur=0.25*q.a_ur(i,ij);
    cl+=dl;cc+=dc;cr+=dr;
    cl+=ul;cc+=uc;cr+=ur;
    dc+=dl;cc+=cl;uc+=ul;
    dc+=dr;cc+=cr;uc+=ur;
    *tp+=dl;tp[1]+=dc;tp[2]+=dr;tp[3]+=cl;tp[4]+=cc;tp[5]+=cr;tp[6]+=ul;tp[7]+=uc;tp[8]=ur;

    // Middle right contribution
    cij=ij+1;ci=i+1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.5*q.a_cl(ci,cij);cc=0.25*q.a_cc(ci,cij);cr=0.5*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    cl+=dl;cc+=dc;cr+=dr;
    cl+=ul;cc+=uc;cr+=ur;
    dl+=dc;dc+=dr;cl+=cc;cc+=cr;ul+=uc;uc+=ur;
    tp[1]+=dl;tp[2]+=dc;tp[4]+=cl;tp[5]+=cc;tp[7]+=ul;tp[8]+=uc;

    // Top left contribution
    cij=ij+m-1;ci=i-1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
    dr+=dc;dc+=dl;cr+=cc;cc+=cl;
    tp[3]+=dc;tp[4]+=dr;tp[6]+=cc;tp[7]+=cr;

    // Top center contribution
    cij=ij+m;
    dl=0.25*q.a_dl(i,cij);dc=0.5*q.a_dc(i,cij);dr=0.25*q.a_dr(i,cij);
    cl=0.125*q.a_cl(i,cij);cc=0.25*q.a_cc(i,cij);cr=0.125*q.a_cr(i,cij);
    ul=0.25*q.a_ul(i,cij);uc=0.5*q.a_uc(i,cij);ur=0.25*q.a_ur(i,cij);
    dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
    dc+=dl;cc+=cl;
    dc+=dr;cc+=cr;
    tp[3]+=dl;tp[4]+=dc;tp[5]+=dr;tp[6]+=cl;tp[7]+=cc;tp[8]+=cr;

    // Top right contribution
    cij=ij+m+1;ci=i+1;
    dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
    cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
    ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
    dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
    dl+=dc;dc+=dr;cl+=cc;cc+=cr;
    tp[4]+=dl;tp[5]+=dc;tp[7]+=cl;tp[8]+=cc;

    // Store reciprocal of central element and update pointer
    tp[9]=1./tp[4];tp+=10;
}

/** Calculates the elements of the matrix at one grid point on the next level
 * down the hierarchy by conjugating with the restriction and interpolation
 * operators, also taking into account modifications due to boundaries.
 * \param[in,out] tp a reference to memory in which to store the matrix
 *                   elements.
 * \param[in] w encoded information about which boundaries need to be accounted
 *              for.
 * \param[in] i the horizontal co-ordinate of the grid point.
 * \param[in] ij the grid point index. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::compute_rat_bdry(M*& tp,unsigned int w,int i,int ij) {
    int ci,cij;
    M dl,dc,dr,cl,cc,cr,ul,uc,ur;

    // Bottom left contribution
    if((w&34)==34) {
        cij=ij-(w&16?m-mn:m)-(w&1?1-m:1);ci=i-(w&1?1-m:1);
        dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
        ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
        cr+=cc;cc+=cl;ur+=uc;uc+=ul;
        *tp=cc;tp[1]=cr;tp[3]=uc;tp[4]=ur;
    } else {*tp=M(0.0);tp[1]=M(0.0);tp[3]=M(0.0);tp[4]=M(0.0);}

    // Bottom center contribution
    if(w&32) {
        cij=ij-(w&16?m-mn:m);
        dl=w&2?0.25*q.a_dl(i,cij):(w&1?0.5*q.a_dl(i,cij):M(0.0));
        dc=0.5*q.a_dc(i,cij);
        dr=w&8?0.25*q.a_dr(i,cij):(w&4?0.5*q.a_dr(i,cij):M(0.0));
        cl=w&2?0.125*q.a_cl(i,cij):(w&1?0.25*q.a_cl(i,cij):M(0.0));
        cc=0.25*q.a_cc(i,cij);
        cr=w&8?0.125*q.a_cr(i,cij):(w&4?0.25*q.a_cr(i,cij):M(0.0));
        ul=w&2?0.25*q.a_ul(i,cij):(w&1?0.5*q.a_ul(i,cij):M(0.0));
        uc=0.5*q.a_uc(i,cij);
        ur=w&8?0.25*q.a_ur(i,cij):(w&4?0.5*q.a_ur(i,cij):M(0.0));
        ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
        if(w&2) {cc+=cl;uc+=ul;}
        if(w&8) {cc+=cr;uc+=ur;}
        *tp+=cl;tp[1]+=cc;tp[2]=cr;tp[3]+=ul;tp[4]+=uc;tp[5]=ur;
    } else {tp[2]=M(0.0);tp[5]=M(0.0);}

    // Bottom right contribution
    if((w&40)==40) {
        cij=ij-(w&16?m-mn:m)+1;ci=i+1;
        dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
        ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        ul+=cl;cl+=dl;uc+=cc;cc+=dc;ur+=cr;cr+=dr;
        cl+=cc;cc+=cr;ul+=uc;uc+=ur;
        tp[1]+=cl;tp[2]+=cc;tp[4]+=ul;tp[5]+=uc;
    }

    // Middle left contribution
    if(w&2) {
        cij=ij-(w&1?1-m:1);ci=i-(w&1?1-m:1);
        if(w&32) {
            dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        } else if(w&16) {
            dl=0.5*q.a_dl(ci,cij);dc=0.25*q.a_dc(ci,cij);dr=0.5*q.a_dr(ci,cij);
        } else {dl=M(0.0);dc=M(0.0);dr=M(0.0);}
        cl=0.5*q.a_cl(ci,cij);cc=0.25*q.a_cc(ci,cij);cr=0.5*q.a_cr(ci,cij);
        if(w&128) {
            ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        } else if(w&64) {
            ul=0.5*q.a_ul(ci,cij);uc=0.25*q.a_uc(ci,cij);ur=0.5*q.a_ur(ci,cij);
        } else {ul=M(0.0);uc=M(0.0);ur=M(0.0);}
        if(w&32) {cl+=dl;cc+=dc;cr+=dr;}
        if(w&128) {cl+=ul;cc+=uc;cr+=ur;}
        dr+=dc;dc+=dl;cr+=cc;cc+=cl;ur+=uc;uc+=ul;
        *tp+=dc;tp[1]+=dr;tp[3]+=cc;tp[4]+=cr;tp[6]=uc;tp[7]=ur;
    } else {tp[6]=M(0.0);tp[7]=M(0.0);}

    // Middle center contribution
    if(w&32) {
        dl=w&2?0.25*q.a_dl(i,ij):(w&1?0.5*q.a_dl(i,ij):M(0.0));
        dc=0.5*q.a_dc(i,ij);
        dr=w&8?0.25*q.a_dr(i,ij):(w&4?0.5*q.a_dr(i,ij):M(0.0));
    } else if(w&16) {
        dl=w&2?0.5*q.a_dl(i,ij):(w&1?q.a_dl(i,ij):M(0.0));
        dc=q.a_dc(i,ij);
        dr=w&8?0.5*q.a_dr(i,ij):(w&4?q.a_dr(i,ij):M(0.0));
    } else {dl=M(0.0);dc=M(0.0);dr=M(0.0);}
    cl=w&2?0.5*q.a_cl(i,ij):(w&1?q.a_cl(i,ij):M(0.0));
    cc=q.a_cc(i,ij);
    cr=w&8?0.5*q.a_cr(i,ij):(w&4?q.a_cr(i,ij):M(0.0));
    if(w&128) {
        ul=w&2?0.25*q.a_ul(i,ij):(w&1?0.5*q.a_ul(i,ij):M(0.0));
        uc=0.5*q.a_uc(i,ij);
        ur=w&8?0.25*q.a_ur(i,ij):(w&4?0.5*q.a_ur(i,ij):M(0.0));
    } else if(w&64) {
        ul=w&2?0.5*q.a_ul(i,ij):(w&1?q.a_ul(i,ij):M(0.0));
        uc=q.a_uc(i,ij);
        ur=w&8?0.5*q.a_ur(i,ij):(w&4?q.a_ur(i,ij):M(0.0));
    } else {ul=M(0.0);uc=M(0.0);ur=M(0.0);}
    if(w&32) {cl+=dl;cc+=dc;cr+=dr;}
    if(w&128) {cl+=ul;cc+=uc;cr+=ur;}
    if(w&2) {dc+=dl;cc+=cl;uc+=ul;}
    if(w&8) {dc+=dr;cc+=cr;uc+=ur;}
    *tp+=dl;tp[1]+=dc;tp[2]+=dr;tp[3]+=cl;tp[4]+=cc;tp[5]+=cr;tp[6]+=ul;tp[7]+=uc;tp[8]=ur;

    // Middle right contribution
    if(w&8) {
        cij=ij+1;ci=i+1;
        if(w&32) {
            dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        } else if(w&16) {
            dl=0.5*q.a_dl(ci,cij);dc=0.25*q.a_dc(ci,cij);dr=0.5*q.a_dr(ci,cij);
        } else {dl=M(0.0);dc=M(0.0);dr=M(0.0);}
        cl=0.5*q.a_cl(ci,cij);cc=0.25*q.a_cc(ci,cij);cr=0.5*q.a_cr(ci,cij);
        if(w&128) {
            ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        } else if(w&64) {
            ul=0.5*q.a_ul(ci,cij);uc=0.25*q.a_uc(ci,cij);ur=0.5*q.a_ur(ci,cij);
        } else {ul=M(0.0);uc=M(0.0);ur=M(0.0);}
        if(w&32) {cl+=dl;cc+=dc;cr+=dr;}
        if(w&128) {cl+=ul;cc+=uc;cr+=ur;}
        dl+=dc;dc+=dr;cl+=cc;cc+=cr;ul+=uc;uc+=ur;
        tp[1]+=dl;tp[2]+=dc;tp[4]+=cl;tp[5]+=cc;tp[7]+=ul;tp[8]+=uc;
    }

    // Top left contribution
    if((w&130)==130) {
        cij=ij+m-(w&1?1-m:1);ci=i-(w&1?1-m:1);
        dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
        ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
        dr+=dc;dc+=dl;cr+=cc;cc+=cl;
        tp[3]+=dc;tp[4]+=dr;tp[6]+=cc;tp[7]+=cr;
    }

    // Top center contribution
    if(w&128) {
        cij=ij+m;
        dl=w&2?0.25*q.a_dl(i,cij):(w&1?0.5*q.a_dl(i,cij):M(0.0));
        dc=0.5*q.a_dc(i,cij);
        dr=w&8?0.25*q.a_dr(i,cij):(w&4?0.5*q.a_dr(i,cij):M(0.0));
        cl=w&2?0.125*q.a_cl(i,cij):(w&1?0.25*q.a_cl(i,cij):M(0.0));
        cc=0.25*q.a_cc(i,cij);
        cr=w&8?0.125*q.a_cr(i,cij):(w&4?0.25*q.a_cr(i,cij):M(0.0));
        ul=w&2?0.25*q.a_ul(i,cij):(w&1?0.5*q.a_ul(i,cij):M(0.0));
        uc=0.5*q.a_uc(i,cij);
        ur=w&8?0.25*q.a_ur(i,cij):(w&4?0.5*q.a_ur(i,cij):M(0.0));
        dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
        if(w&2) {dc+=dl;cc+=cl;}
        if(w&8) {dc+=dr;cc+=cr;}
        tp[3]+=dl;tp[4]+=dc;tp[5]+=dr;tp[6]+=cl;tp[7]+=cc;tp[8]+=cr;
    }

    // Top right contribution
    if((w&136)==136) {
        cij=ij+m+1;ci=i+1;
        dl=0.25*q.a_dl(ci,cij);dc=0.125*q.a_dc(ci,cij);dr=0.25*q.a_dr(ci,cij);
        cl=0.125*q.a_cl(ci,cij);cc=0.0625*q.a_cc(ci,cij);cr=0.125*q.a_cr(ci,cij);
        ul=0.25*q.a_ul(ci,cij);uc=0.125*q.a_uc(ci,cij);ur=0.25*q.a_ur(ci,cij);
        dl+=cl;cl+=ul;dc+=cc;cc+=uc;dr+=cr;cr+=ur;
        dl+=dc;dc+=dr;cl+=cc;cc+=cr;
        tp[4]+=dl;tp[5]+=dc;tp[7]+=cl;tp[8]+=cc;
    }

    // Store reciprocal of central element and update pointer
    tp[9]=1./tp[4];tp+=10;
}

/** Copies the source array into the solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::b_to_z() {
#pragma omp parallel for num_threads(nthr_simp)
    for(int j=0;j<mn;j+=m) {
        for(V *bp=b+j,*zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=*(bp++);
    }
}

/** Clears the solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::clear_z() {
#pragma omp parallel for num_threads(nthr_simp)
    for(int j=0;j<mn;j+=m) {
        for(V *zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=V(0.);
    }
}

/** Carries out a Jacobi iteration assuming that the current solution is zero.
 * This is a quick operation and may provide a better initial guess for the
 * Gauss--Seidel sweeps than setting the solution to zero. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::zero_jacobi() {
#pragma omp parallel for num_threads(nthr_simp)
    for(int j=0;j<mn;j+=m) {
        for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]);
    }
}

/** Saves a field to a file in the gnuplot matrix binary format. This is mainly
 * used for debugging and diagnostic purposes.
 * \param[in] filename the name of the file to save to.
 * \param[in] ff a pointer to the field to save. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::output(const char *filename,V *ff,double ax,double dx,double ay,double dy) {

    // Assemble the output filename and open the output file
    float *buf=new float[m+1];
    FILE *outf=fopen(filename,"wb");
    if(outf==NULL) {
        fputs("File output error\n",stderr);
        exit(TGMGPP_ERROR);
    }

    // Output the first line of the file
    int i,j;
    *buf=m;
    float *bp=buf+1,*be=buf+m+1;
    for(i=0;i<m;i++) *(bp++)=ax+i*dx;
    fwrite(buf,sizeof(float),m+1,outf);

    // Output the field values to the file
    V *fp=ff;
    for(j=0;j<n;j++) {
        *buf=ay+j*dy;bp=buf+1;
        while(bp<be) *(bp++)=float_output(*(fp++));
        fwrite(buf,sizeof(float),m+1,outf);
    }

    // Close the file
    fclose(outf);
    delete [] buf;
}

/** Outputs the matrix residual to a file.
 * \param[in] filename the name of the file to save to. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::output_res(const char* filename,double ax,double dx,double ay,double dy) {
    V *r=new V[mn];
#pragma omp parallel for num_threads(nthr_r)
    for(int j=0;j<mn;j+=m) {
        for(int i=0,ij=j;i<m;i++,ij++) r[ij]=res(i,ij);
    }
    output(filename,r,ax,dx,ax,dy);
    delete [] r;
}
