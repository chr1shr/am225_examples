#include <limits>

#include "tgmg.hh"

/** The multigrid constructor determines constants from the multigrid setup
 * class, and creates the hierarchy of coarser grids.
 * \param[in] q_ the multigrid setup class to use.
 * \param[in] b_ a pointer to the source array.
 * \param[in] z_ a pointer to the solution array. */
template<class S,class V,class M>
tgmg<S,V,M>::tgmg(S &q_,V* b_,V* z_) : tgmg_base<S,V,M>(q_,q_.m,q_.n,q_.x_prd,q_.y_prd,b_,z_), ml(0),
	verbose(1), conv_rate(std::numeric_limits<double>::max()) {
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
	printf("Top grid level: (%d,%d) [%s,%s] {%d}\n",
	       m,n,m&1?"odd":"even",n&1?"odd":"even",num_t);
	for(int l=0;l<ml;l++) {
		int am=mg[l]->sm,an=mg[l]->sn;
		printf("Grid level %2d : (%d,%d) [%s,%s] {%d}\n",
		       ml,am,an,am&1?"odd":"even",an&1?"odd":"even",mg[l]->num_t);
	}
}

/** Calculates the sum of squares of residuals.
 * \return The sum. */
template<class S,class V,class M>
double tgmg_base<S,V,M>::mds() {
	double c=0;
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		double cl=0;
		for(int i=0,ij=j;i<m;i++,ij++) cl+=mod_sq(b[ij]-q.mul_a(i,ij)-q.a_cc(i,ij)*z[ij]);
		#pragma omp atomic
		c+=cl;
	}
	return c;
}

/** Carries out a Jacobi iteration, by creating a temporary array to store the
 * improved solution, and then copying it over the main solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::jacobi() {
	V *zn=new V[mn];
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		for(int i=0,ij=j;i<m;i++,ij++) zn[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
	}
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		for(V *znp=zn+j,*zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=*(znp++);
	}
	delete [] zn;
}

/** Carries out a Gauss-Seidel sweep. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::gauss_seidel() {
	if(q.gs_mode<2) {

		// Regular Gauss--Seidel
#pragma omp parallel for num_threads(num_t)
		for(int j=0;j<mn;j+=2*m) {
			for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
		}
#pragma omp parallel for num_threads(num_t)
		for(int j=m;j<mn;j+=2*m) {
			for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
		}

		// Additional step for symmetric Gauss--Seidel
		if(q.gs_mode==1) {
#pragma omp parallel for num_threads(num_t)
			for(int j=2*m-1;j<mn;j+=2*m) {
				for(int i=m-1,ij=j;i>=0;i--,ij--) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
			}
#pragma omp parallel for num_threads(num_t)
			for(int j=m-1;j<mn;j+=2*m) {
				for(int i=m-1,ij=j;i>=0;i--,ij--) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
			}
		}
	} else {

		// Crazy Gauss-Seidel where threads don't care about about if
		// they're accessing new or old data. I'm not sure if this is
		// legit or not. If writing V entries is not an atomic
		// operation then bad things could happen. But most of the time
		// it's probably fine.
#pragma omp parallel for num_threads(num_t)
		for(int j=0;j<mn;j+=m) {
			for(int i=0,ij=j;i<m;i++,ij++) z[ij]=q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij));
		}
	}
}

/** Carries out a successive over-relation (SOR) sweep.
 * \param[in] omega the factor by which to relax by, usually between 0 and 2.
 */
template<class S,class V,class M>
void tgmg_base<S,V,M>::sor(double omega) {
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=2*m) {
		for(int i=0,ij=j;i<m;i++,ij++) z[ij]+=omega*(q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij))-z[ij]);
	}
#pragma omp parallel for num_threads(num_t)
	for(int j=m;j<mn;j+=2*m) {
		for(int i=0,ij=j;i<m;i++,ij++) z[ij]+=omega*(q.inv_cc(i,ij,b[ij]-q.mul_a(i,ij))-z[ij]);
	}
}

/** Solves the linear system using one of the smoothing techniques. It performs
 * batches of smoothing cycles and checks after each to see if the specified
 * tolerance is reached, after which it terminates.
 * \param[in] type the type of smoothing to perform.
 * \param[in] per_loop the number of smoothing cycles to perform in each batch.
 * \param[in] omega the over-relaxation parameter, only used with SOR smoothing.
 * \param[in] (cyc_down,cyc_up,cyc_bottom,cyc_top)
 *            V-cycle configuration parameters.
 * \return True if the tolerance was reached before the maximum number of
 * iterations, false otherwise. */
template<class S,class V,class M>
bool tgmg<S,V,M>::solve(int type,int per_loop,int max_loops,double omega,int cyc_down,int cyc_up,int cyc_bottom,int cyc_top) {
	double iacc=1,acc=0;int n=0;
	if(verbose>=2) iacc=l2_error();
	do {
		if(++n==max_loops) {
			bail_message(n*per_loop,iacc,acc);
			return false;
		}
		acc=iters_and_error(type,per_loop,omega,cyc_down,cyc_up,cyc_bottom,cyc_top);
		if(verbose==3) iter_message(n*per_loop,acc);
	} while(acc>q.acc);
	status_message(n*per_loop,iacc,acc);
	return !std::isnan(acc);
}

/** Solves the linear system using one of the smoothing techniques and the adaptive
 * approach for choosing iterations.
 * \param[in] type the type of smoothing to perform.
 * \param[in] tp a class for predicting the number of smoothing steps needed.
 * \param[in] omega the over-relaxation parameter, only used with SOR smoothing.
 * \param[in] (cyc_down,cyc_up,cyc_bottom,cyc_top)
 *            V-cycle configuration parameters.
 * \return True if the tolerance was reached before the maximum number of
 * iterations, false otherwise. */
template<class S,class V,class M>
bool tgmg<S,V,M>::solve(int type,tgmg_predict &tp,double omega,int cyc_down,int cyc_up,int cyc_bottom,int cyc_top) {
	int n=tp.lim/tp.mult;

	// Perform the predicted number of smoothing iterations
	double iacc=verbose>=2?l2_error():1,
	       acc=iters_and_error(type,n,omega,cyc_down,cyc_up,cyc_bottom,cyc_top);
	if(verbose==3) iter_message(n,acc);

	// Test the L2 error against the given threshold. If it's lower, then
	// try to decrease the smoothing iterations for the solve.
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
			acc=iters_and_error(type,k,omega,cyc_down,cyc_up,cyc_bottom,cyc_top);
			if(verbose==3) iter_message(n,acc);
		} while(acc>=q.acc);
	}

	// Record the number of iterations and perform any extra smoothing
	// steps
	tp.add_iters(n);
	status_message(n,iacc,acc);
	iters(type,tp.extra_iters,omega,cyc_down,cyc_up,cyc_bottom,cyc_top);
	return !std::isnan(acc);
}

/** Carries out a V-cycle.
 * \param[in] cyc_down the number of Gauss--Seidel sweeps to apply on the way
 *		       down the grid hierarchy.
 * \param[in] cyc_up the number of Gauss--Seidel sweeps to apply on the way up
 *		     the grid hierarchy, not including the top level.
 * \param[in] cyc_bottom the number of Gauss--Seidel sweeps to apply on the
 *			 bottom level of the grid hierarchy.
 * \param[in] cyc_top the number of Gauss--Seidel sweeps to apply on the top
 *		      level of the grid hierarchy. */
template<class S,class V,class M>
void tgmg<S,V,M>::v_cycle(int cyc_down,int cyc_up,int cyc_bottom,int cyc_top) {
	int i,j;

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

		// Propagate the solution up the hierarchy, smoothing at each step
		for(i=ml-1;i>=0;i--) {
			if(i!=ml-1) {
				for(j=0;j<cyc_up;j++) mg[i]->gauss_seidel();
			}
			mg[i]->apply_t();
		}
	}

	// Apply smoothing sweeps on the top level
	for(j=0;j<cyc_top;j++) gauss_seidel();
}

/** Prints the calculated matrix entries. This function is mainly used for
 * diagnostic purposes and debugging. */
template<class V,class M>
void tgmg_level<V,M>::print_rat() {
	int i,j,k;
	M* sp=s;
	for(j=0;j<n;j++) {
		for(i=0;i<m;i++,sp+=10) {
			printf("%d %d [",i,j);
			for(k=0;k<8;k++) printf("%g,",sp[k]);
			printf("%g] {%g}\n",sp[8],sp[9]);
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
	if(ij<m) {
		if (i==0) {
			return tsub_contrib(0,0);
		} else if (i==m-1) {
			return tsub_contrib(2,0);
		} else {
			return tsub_contrib(1,0);
		}
	} else if (ij>=mn-m) {
		if (i==0) {
			return tsub_contrib(0,2);
		} else if (i==m-1) {
			return tsub_contrib(2,2);
		} else {
			return tsub_contrib(1,2);
		}
	} else {
		if (i==0) {
			return tsub_contrib(0,1);
		} else if (i==m-1) {
			return tsub_contrib(2,1);
		} else {
			return tsub_contrib(1,1);
		}
	}
}

/** Applies the restriction operation, computing the residual at the current
 * level and projecting it the level below. */
template<class V,class M>
void tgmg_level<V,M>::apply_t() {

	// Deal with the bulk of the grid
#pragma omp parallel for num_threads(num_t)
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

	// Handle the first bulk pass where the residuals are set into the
	// child array
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<n-2;j+=4) r_double_line_set(c+sm*(j>>1),j*m);

	// Handle the boundary lines, taking into account periodicity if
	// necessary
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
#pragma omp parallel for num_threads(num_t)
	for(int j=2;j<n-2;j+=4) r_double_line_add(c+sm*(j>>1),j*m);
}

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_double_line_set(V* cp,int ij);

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_double_line_add(V* cp,int ij);

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_line_set(V* cp,int ij);

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_line_add(V* cp,int ij);

template<class S,class V,class M>
void tgmg_base<S,V,M>::r_periodic_line_add(V* cp,int ij);

/** Calculates the matrix entries on the child grid by conjugating the matrix
 * entries on this grid by the restriction and interpolation operators. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::rat() {
	if(y_prd) {
		if((n&1)==0) {
			rat_line(t,176,0); //IR
			if(n>2) {
#pragma omp parallel for num_threads(num_t)
				for(int j=2;j<n;j+=2) rat_line(t+5*sm*j,160,j*m);
			} else rat_collapse_y();
		} else if(n==1) rat_line(t,0,0); //HH
		else {
			rat_line(t,144,0); //DR
#pragma omp parallel for num_threads(num_t)
			for(int j=2;j<n-2;j+=2) rat_line(t+5*sm*j,160,j*m); //RR
			rat_line(t+10*sm*(sn-1),96,(n-1)*m); //RD
		}
	} else {
		if(n==2) {
			rat_line(t,64,0); //HD
			rat_line(t+10*sm,16,1); //DH
		} else {
			rat_line(t,128,0); //HR
#pragma omp parallel for num_threads(num_t)
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
void tgmg_base<S,V,M>::compute_rat(M*& tp,int i,int ij);

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
void tgmg_base<S,V,M>::compute_rat_bdry(M*& tp,unsigned int w,int i,int ij);

/** Copies the source array into the solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::b_to_z() {
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		for(V *bp=b+j,*zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=*(bp++);
	}
}

/** Clears the solution array. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::clear_z() {
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		for(V *zp=z+j,*ze=zp+m;zp<ze;) *(zp++)=0;
	}
}

/** Carries out a Jacobi iteration assuming that the current solution is zero.
 * This is a quick operation and provides a better initial guess for the
 * Gauss--Seidel sweeps than setting the solution to zero. */
template<class S,class V,class M>
void tgmg_base<S,V,M>::zero_jacobi() {
#pragma omp parallel for num_threads(num_t)
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
#pragma omp parallel for num_threads(num_t)
	for(int j=0;j<mn;j+=m) {
		for(int i=0,ij=j;i<m;i++,ij++) r[ij]=res(i,ij);
	}
	output(filename,r,ax,dx,ax,dy);
	delete [] r;
}
