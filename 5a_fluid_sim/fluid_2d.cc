#include <cstring>
#include <limits>

#include "common.hh"
#include "fluid_2d.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] (x_prd_,y_prd_) the periodicity in the x and y directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] visc_ the fluid viscosity.
 * \param[in] rho_ the fluid density.
 * \param[in] filename_ the filename of the output directory. */
fluid_2d::fluid_2d(const int m_,const int n_,const bool x_prd_,
        const bool y_prd_,const double ax_,const double bx_,
        const double ay_,const double by_,const double visc_,
        const double rho_,unsigned int fflags_,const char *filename_)
    : m(m_), n(n_), mn(m_*n_), m_fem(x_prd_?m:m+1), n_fem(y_prd_?n:n+1),
    ml(m+4), ntrace(0), x_prd(x_prd_), y_prd(y_prd_), ax(ax_), ay(ay_),
    bx(bx_), by(by_), dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx), ysp(1/dy),
    xxsp(xsp*xsp), yysp(ysp*ysp), visc(visc_), rho(rho_), rhoinv(1/rho),
    filename(filename_), fbase(new field[ml*(n+4)]), fm(fbase+2*ml+2),
    src(new double[m_fem*n_fem]), tm(NULL), time(0.), f_num(0),
    fflags(fflags_), ms_fem(*this), buf(new float[m>123?m+5:128]) {}

/** The class destructor frees the dynamically allocated memory. */
fluid_2d::~fluid_2d() {
    if(ntrace>0) delete [] tm;
    delete [] buf;
    delete [] src;
    delete [] fbase;
}

/** Initializes the simulation, setting up the tracers and simulation fields,
 * and choosing the timestep.
 * \param[in] ntrace_ the number of tracers.
 * \param[in] dt_pad_ the padding factor for the timestep, which should be
 *                    smaller than 1.
 * \param[in] max_spd a maximum fluid speed from which to estimate the
 *                    advection timestep restriction. If a negative value is
 *                    supplied, then the advection CFL condition is explicitly
 *                    calculated. */
void fluid_2d::initialize(int ntrace_,double dt_pad,double max_spd) {

    // Set up the tracers (if any) and initialize the simulation fields
    if((ntrace=ntrace_)>0) init_tracers();
    init_fields();

    // Compute the timestep, based on the restrictions from the advection
    // and velocity, plus a padding factor
    choose_dt(dt_pad,max_spd<=0?advection_dt():(dx>dy?dx:dy)/max_spd);
}

/** Computes the maximum timestep that can resolve the fluid advection, based
 * on the CFL condition.
 * \return The maximum timestep. */
double fluid_2d::advection_dt() {
    double adv_dt=0;
#pragma omp parallel for reduction(max:adv_dt)
    for(int j=0;j<n;j++) {
        double t;
        for(field *fp=fm+ml*j,*fe=fp+m;fp<fe;fp++) {
            t=fp->cfl(xsp,ysp);
            if(t>adv_dt) adv_dt=t;
        }
    }
    return adv_dt==0?std::numeric_limits<double>::max():1./adv_dt;
}

/** Chooses the timestep based on the limits from advection and viscosity.
 * \param[in] dt_pad the padding factor for the timestep for the physical
 *                   terms, which should be smaller than 1.
 * \param[in] adv_dt the maximum timestep to resolve the fluid advection.
 * \param[in] verbose whether to print out messages to the screen. */
void fluid_2d::choose_dt(double dt_pad,double adv_dt,bool verbose) {

    // Calculate the viscous timestep restriction
    double vis_dt=0.5*rho/(visc*(xxsp+yysp));
    int ca;

    // Choose the minimum of the two timestep restrictions
    if(adv_dt<vis_dt) {ca=0;dt_reg=adv_dt;}
    else {ca=1;dt_reg=vis_dt;}
    dt_reg*=dt_pad;

    // Print information if requested
    if(verbose) {
        const char mno[]="", myes[]=" <-- use this";
        printf("# Advection dt       : %g%s\n"
               "# Viscous dt         : %g%s\n"
               "# Padding factor     : %g\n"
               "# Minimum dt         : %g\n",
               adv_dt,ca==0?myes:mno,vis_dt,ca==1?myes:mno,dt_pad,dt_reg);
    }
}

/** Initializes the simulation fields. */
void fluid_2d::init_fields() {

    // Loop over the primary grid and set the velocity and pressure
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay+dy*(j+0.5),yy=y+0.5;
        field *fp=fm+ml*j;
        for(int i=0;i<m;i++) {
            double x=ax+dx*(i+0.5),xx=x+0.5;
            fp->u=4*exp(-20*(x*x+y*y));
            fp->v=exp(-20*(x*x+y*y))-5*exp(-30*(xx*xx+yy*yy));
            (fp++)->p=0;
        }
        fp->p=0;
    }

    // Set the final line of the cell-cornered pressure field
    field *fp=fm+ml*n,*fe=fp+m+1;
    while(fp<fe) (fp++)->p=0;

    // Now that the primary grid points are set up, initialize the ghost
    // points according to the boundary conditions
    set_boundaries();
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void fluid_2d::solve(double duration,int frames) {
    double t0,t1,t2,adt;
    int l=timestep_select(duration/frames,adt);

    // Save header file, output the initial fields and record the initial wall
    // clock time
    save_header(duration,frames);
    if(f_num==0) write_files(0), puts("# Output frame 0");
    t0=wtime();

    // Loop over the output frames
    for(int k=1;k<=frames;k++) {

        // Perform the simulation steps
        for(int j=0;j<l;j++) step_forward(adt);

        // Output the fields
        t1=wtime();
        write_files(k+f_num);

        // Print diagnostic information
        t2=wtime();
        printf("# Output frame %d [%d, %.8g s, %.8g s] {FEM %.2f}\n",
               k,l,t1-t0,t2-t1,ms_fem.tp.avg_iters());
        t0=t2;
    }
    f_num+=frames;
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void fluid_2d::step_forward(double dt) {
    int j;
    double hx=0.5*dt*xsp,hy=0.5*dt*ysp,hxx=rhoinv*visc*xxsp*dt,
           hyy=rhoinv*visc*yysp*dt;

    // Perform an explicit Euler step of the tracer positions using the
    // bilinear interpolation of the velocity field
    update_tracers(dt);

#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {
        field *fp=fm+(ml*j+i),&f=*fp;

        // Compute the second derivatives that are needed to evaluate the
        // viscous stresses
        double ux,vx,uy,vy,&uc=f.u,&vc=f.v,
               uyy=hyy*(fp[-ml].u-2*uc+fp[ml].u),
               vyy=hyy*(fp[-ml].v-2*vc+fp[ml].v),
               uxx=hxx*(fp[-1].u-2*uc+fp[1].u),
               vxx=hxx*(fp[-1].v-2*vc+fp[1].v);

        // Compute advective terms using the second-order ENO scheme
        uc>0?vel_eno2(ux,vx,hx,fp[1],f,fp[-1],fp[-2])
            :vel_eno2(ux,vx,-hx,fp[-1],f,fp[1],fp[2]);
        vc>0?vel_eno2(uy,vy,hy,fp[ml],f,fp[-ml],fp[-2*ml])
            :vel_eno2(uy,vy,-hy,fp[-ml],f,fp[ml],fp[2*ml]);

        // Compute the intermediate velocity using advection and viscosity.
        // Note that the terms ux, uyy, etc. are already scaled by the correct
        // constants.
        f.us=f.u-uc*ux-vc*uy+uxx+uyy;
        f.vs=f.v-uc*vx-vc*vy+vxx+vyy;
    }

    // Calculate the source term for the finite-element projection, doing
    // some preliminary copying to simplify periodicity calculations
    fem_source_term_conditions();
    double sx=0.5*dx/dt,sy=0.5*dy/dt;
#pragma omp parallel for
    for(j=0;j<n_fem;j++) {
        double *srp=src+j*m_fem;
        for(field *fp=fm+j*ml,*fe=fp+m_fem;fp<fe;fp++)
            *(srp++)=sx*(fp[-ml-1].us+fp[-1].us-fp[-ml].us-fp->us)
                    +sy*(fp[-ml-1].vs-fp[-1].vs+fp[-ml].vs-fp->vs);
    }

    // Solve the finite-element problem, and copy the pressure back into
    // the main data structure, subtracting off the mean, and taking into
    // account the boundary conditions
    ms_fem.solve_v_cycle();
    copy_pressure();

    // Update u and v based on us, vs, and the computed pressure
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        while(fp<fe) {
            fp->u=fp->us-0.5*dt*rhoinv*xsp*(fp[ml+1].p+fp[1].p-fp[ml].p-fp->p);
            fp->v=fp->vs-0.5*dt*rhoinv*ysp*(fp[ml+1].p-fp[1].p+fp[ml].p-fp->p);
            fp++;
        }
    }

    // Reset the ghost points according to the boundary conditions
    set_boundaries();

    // Increment time at end of step
    time+=dt;
}

/** Copies the pressure back into the main data structure, subtracting off the
 * mean, and taking into account the boundary conditions. */
void fluid_2d::copy_pressure() {
    double pavg=average_pressure(),*sfem=ms_fem.z;

#pragma omp parallel for
    for(int j=0;j<n+1;j++) {

        // Set a pointer to the row to copy. If the domain is
        // y-periodic and this is the last line, then
        double *sop=y_prd&&j==n?sfem:sfem+j*m_fem,*soe=sop+m;
        field *fp=fm+j*ml;
        while(sop<soe) (fp++)->p=*(sop++)-pavg;
        fp->p=x_prd?fp[-m].p:*sop-pavg;
    }
}

/** Computes the average pressure that has been computed using the
 * finite-element method.
 * \return The pressure. */
double fluid_2d::average_pressure() {
    double pavg=0,*sfem=ms_fem.z;

#pragma omp parallel for reduction(*:pavg)
    for(int j=0;j<n_fem;j++) {
        double *sop=sfem+j*m_fem,*soe=sop+(m_fem-1);
        double prow=*(sop++)*(x_prd?1:0.5);
        while(sop<soe) prow+=*(sop++);
        prow+=*sop*(x_prd?1:0.5);
        if(!y_prd&&(j==0||j==n)) prow*=0.5;
        pavg+=prow;
    }
    return pavg*(1./mn);
}

/** Calculates one-sided derivatives of the velocity field using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[out] (ud,vd) the computed ENO2 derivatives.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with. */
inline void fluid_2d::vel_eno2(double &ud,double &vd,double hs,field &f0,field &f1,field &f2,field &f3) {
    ud=hs*eno2(f0.u,f1.u,f2.u,f3.u);
    vd=hs*eno2(f0.v,f1.v,f2.v,f3.v);
}

/** Calculates the ENO derivative using a sequence of values at four
 * gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double fluid_2d::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void fluid_2d::set_boundaries() {

    // Set left and right ghost values
    if(x_prd) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].prd_bc(fp[m-2]);
            fp[-1].prd_bc(fp[m-1]);
            fp[m].prd_bc(*fp);
            fp[m+1].prd_bc(fp[1]);
        }
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].no_slip(fp[1]);
            fp[-1].no_slip(*fp);
            fp[m].no_slip(fp[m-1]);
            fp[m+1].no_slip(fp[m-2]);
        }
    }

    // Set top and bottom ghost values
    const int tl=2*ml,g=n*ml;
    if(y_prd) {
        for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) {
            fp[-tl].prd_bc(fp[g-tl]);
            fp[-ml].prd_bc(fp[g-ml]);
            fp[g].prd_bc(*fp);
            fp[g+ml].prd_bc(fp[ml]);
        }
    } else {
        for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) {
            fp[-tl].no_slip(fp[ml]);
            fp[-ml].no_slip(*fp);
            fp[g].no_slip(fp[g-ml]);
            fp[g+ml].no_slip(fp[g-tl]);
        }
    }
}

/** Sets boundary conditions for the FEM source term computation, taking into
 * account periodicity */
void fluid_2d::fem_source_term_conditions() {

    // Set left and right ghost values
    int xl;
    if(x_prd) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].us=fp[m-1].us;fp[-1].vs=fp[m-1].vs;
        }
        xl=m;
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].us=fp[-1].vs=0;
            fp[m].us=fp[m].vs=0;
        }
        xl=m+1;
    }

    // Set top and bottom ghost values
    const int g=n*ml;
    if(y_prd) {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].us=fp[g-ml].us;fp[-ml].vs=fp[g-ml].vs;
        }
    } else {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].us=fp[-ml].vs=0;
            fp[g].us=fp[g].vs=0;
        }
    }
}

/** Sets up the fluid tracers by initializing them at random positions. */
void fluid_2d::init_tracers() {
    tm=new double[ntrace<<1];
    for(double *tp=tm;tp<tm+(ntrace<<1);) {

        // Create a random position vector within the simulation region
        *(tp++)=ax+(bx-ax)/RAND_MAX*double(rand());
        *(tp++)=ay+(by-ay)/RAND_MAX*double(rand());
    }
}

/** Moves the tracers according to the bilinear interpolation of the fluid
* velocity.
* \param[in] dt the timestep to use. */
void fluid_2d::update_tracers(double dt) {
    int i,j;
    double x,y;
    for(double *tp=tm,*te=tm+(ntrace<<1);tp<te;tp+=2) {

        // Find which grid cell the tracer is in
        x=(*tp-ax)*xsp+0.5;y=(tp[1]-ay)*ysp+0.5;
        i=int(x)-1;
        j=int(y)-1;
        if(i<-1) i=-1;else if(i>=m) i=m-1;
        if(j<-1) j=-1;else if(j>=n) i=n-1;

        // Compute the tracer's fractional position with the grid cell
        x-=i;y-=j;

        // Compute tracer's new position
        field *fp=fm+(i+ml*j);
        *tp+=dt*((1-y)*(fp->u*(1-x)+fp[1].u*x)+y*(fp[ml].u*(1-x)+fp[ml+1].u*x));
        tp[1]+=dt*((1-y)*(fp->v*(1-x)+fp[1].v*x)+y*(fp[ml].v*(1-x)+fp[ml+1].v*x));
        remap_tracer(*tp,tp[1]);
    }
}

/** Remap a tracer according to periodic boundary conditions, if they
 * are being used.
 * \param[in,out] (xx,yy) the tracer position, which is updated in situ. */
inline void fluid_2d::remap_tracer(double &xx,double &yy) {
    const double xfac=1./(bx-ax),yfac=1./(by-ay);
    if(x_prd) xx-=(bx-ax)*floor((xx-ax)*xfac);
    if(y_prd) yy-=(by-ay)*floor((yy-ay)*yfac);
}

/** Outputs the tracer positions in a binary format that can be read by
 * Gnuplot. */
void fluid_2d::output_tracers(const char *prefix,const int sn) {
    if(ntrace==0) return;

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the tracer positions in batches of 128 floats
    int j,tbatch=(ntrace<<1)/128,tres=(ntrace<<1)%128;
    float *fp,*fe=buf+128;
    double *tp=tm,*te=tm+(ntrace<<1);
    for(j=0;j<tbatch;j++) {
        fp=buf;
        while(fp<fe) *(fp++)=*(tp++);
        fwrite(buf,sizeof(float),128,outf);
    }

    // Output the remaining tracer positions, if any
    if(tres>0) {
        fp=buf;
        do {*(fp++)=*(tp++);} while(tp<te);
        fwrite(buf,sizeof(float),tres,outf);
    }

    // Close the file
    fclose(outf);
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void fluid_2d::write_files(int k) {
    if(fflags&1) output("u",0,k);
    if(fflags&2) output("v",1,k);
    if(fflags&4) output("p",2,k);
    output_tracers("trace",k);
}

/** Saves the header file.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void fluid_2d::save_header(double duration, int frames) {
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/header",filename);
    FILE *outf=safe_fopen(bufc,f_num==0?"w":"a");
    fprintf(outf,"%g %g %d\n",time,time+duration,frames);
    fclose(outf);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] ghost whether to output the ghost regions or not. */
void fluid_2d::output(const char *prefix,const int mode,const int sn,const bool ghost) {

    // Determine whether to output a cell-centered field or not
    bool cen=mode>=0&&mode<=1;
    int lx=ghost?ml:(cen?m:m+1),ly=ghost?n+4:(cen?n:n+1);
    double disp=(cen?0.5:0)-(ghost?2:0);

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    int i,j;
    float *bp=buf+1,*be=bp+lx;
    *buf=lx;
    for(i=0;i<lx;i++) *(bp++)=ax+(i+disp)*dx;
    fwrite(buf,sizeof(float),lx+1,outf);

    // Output the field values to the file
    field *fr=ghost?fbase:fm;
    for(j=0;j<ly;j++,fr+=ml) {
        field *fp=fr;
        *buf=ay+(j+disp)*dy;bp=buf+1;
        switch(mode) {
            case 0: while(bp<be) *(bp++)=(fp++)->u;break;
            case 1: while(bp<be) *(bp++)=(fp++)->v;break;
            case 2: while(bp<be) *(bp++)=(fp++)->p;break;
        }
        fwrite(buf,sizeof(float),lx+1,outf);
    }

    // Close the file
    fclose(outf);
}
