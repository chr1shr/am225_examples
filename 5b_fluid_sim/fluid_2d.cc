#include <cstring>

#include "common.hh"
#include "fluid_2d.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *              vertical directions.
 * \param[in] (x_prd_,y_prd_) the periodicities in the x and y directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] visc_ the viscosity.
 * \param[in] rho_ the fluid density.
 * \param[in] tmult_ a multiplier to apply to the default timestep size.
 * \param[in] ntrace_ the number of tracers to use.
 * \param[in] filename_ the filename of the output directory. */
fluid_2d::fluid_2d(const int m_,const int n_,const bool x_prd_,const bool y_prd_,
        const double ax_,const double bx_,
        const double ay_,const double by_,
        const double visc_,const double rho_,const double tmult_,
        const int ntrace_,const bool implicit_visc,const char *filename_)
    : m(m_), n(n_), mn(m_*n_), m_fem(x_prd_?m:m+1), n_fem(y_prd_?n:n+1),
    ml(m+4), ntrace(ntrace_), x_prd(x_prd_), y_prd(y_prd_),
    ax(ax_), ay(ay_), bx(bx_), by(by_), dx((bx_-ax_)/m_), dy((by_-ay_)/n_),
    xsp(1/dx), ysp(1/dy), xxsp(xsp*xsp), yysp(ysp*ysp), visc(visc_), rho(rho_),
    rhoinv(1/rho), tmult(tmult_), filename(filename_),
    fbase(new field[ml*(n+4)]), fm(fbase+2*ml+2),
    tm(new double[ntrace<<1]), src(new double[m_fem*n_fem]),
    s_fem(new double[m_fem*n_fem]),
    m_fem(*this), mg1(m_fem,src,s_fem),
    buf(new float[m>123?m+5:128]) {

    // Initialize the multigrid hierarchies
    mg1.setup();mg1.clear_z();

    // Initialize the tracers
    if(ntrace>0) init_tracers();
}

/** The class destructor frees the dynamically allocated memory. */
fluid_2d::~fluid_2d() {
    delete [] buf;
    delete [] sfem;delete [] src;
    delete [] tm;delete [] fbase;
}

/** Initializes the simulation fields. */
void fluid_2d::init_fields() {

#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay+dy*(j+0.5);
        double xx,yy;
        field *fp=fm+ml*j;
        for(int i=0;i<m;i++) {
            double x=ax+dx*(i+0.5);
            fp->u = sin(M_PI*x)*cos(M_PI*y);
            fp->v = -cos(M_PI*x)*sin(M_PI*y);
            fp->p=0;
            fp++;
        }
        double x=ax+dx*(m+0.5);
    }

    // Now that the primary grid points are set up, initialize the ghost
    // points according to the boundary conditions
    set_boundaries();
}

/** Carries out the simulation for a specified duration, periodically saving
 * the output.
 * \param[in] duration the simulation duration
 * \param[in] frames the number of frames to save. */
void fluid_2d::solve(double t_start,double t_end,int frames) {
    time=t_start;
    double time_interval=(t_end-t_start)/frames,t0,t1,t2;
    //const double dt=dx*dx*tmult;
    const double dt=dx*tmult;
    int l=int(time_interval/dt)+1;
    double adt=time_interval/l;

    // Save header file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/header",filename);
    FILE *outf=safe_fopen(bufc,"w");
    fprintf(outf,"%g %g %d\n",t_start,t_end,frames);
    fclose(outf);

    // Output the initial fields and record initial time
    set_boundaries();
    write_files(0);
    puts("# Output frame 0");
    t0=wtime();

    for(int k=1;k<=frames;k++) {

        for(int qq=0;qq<l;qq++) step_forward(adt);
        //time+=time_interval; //BUG

        // Output the fields
        t1=wtime();
        write_files(k);

        // Print diagnostic information
        t2=wtime();
        printf("# Output frame %d [%d, %.8g s, %.8g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void fluid_2d::step_forward(double dt) {
    int j;

    // Update the simulation time, and move the tracers
    update_tracers(dt);

#pragma omp parallel for
    for(j=0;j<n;j++) for(i=0;i<m;i++) {

    }

#pragma omp parallel for
    for(field *fp=fm;fp<fm+mn;fp++) fp->update();

    // compute intermediate velocities at the next time step
#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {

        // Create references to the fields in the neighboring gridpoints
        field *fp=fm+(ml*j+i),&f=*fp,&fr=fp[1],&fu=fp[ml];

                if (disable_nonlinear) {
                  // Disable nonlinear terms
          f.c0=f.u+dt*f.c0;
          f.c1=f.v+dt*f.c1;
                }
                else {
          f.c0=f.u+dt*(-0.5*(xsp*(fr.ul+f.ul)*(fr.ul-f.ul)+ysp*(fu.vd+f.vd)*(fu.ud-f.ud))+f.c0);
          f.c1=f.v+dt*(-0.5*(xsp*(fr.ul+f.ul)*(fr.vl-f.vl)+ysp*(fu.vd+f.vd)*(fu.vd-f.vd))+f.c1);
                }
    }

    // Calculate the source term for the finite-element projection, doing
    // some preliminary copying to simplify periodicity calculations
        fem_source_term_conditions();
        double hx=0.5*dx/dt,hy=0.5*dy/dt;
    #pragma omp parallel for
        for(j=0;j<n_fem;j++) {
            double *srp=src+j*m_fem;
            for(field *fp=fm+j*ml,*fe=fp+m_fem;fp<fe;fp++) {
                *(srp++)=hx*(fp[-ml-1].c0+fp[-1].c0-fp[-ml].c0-fp->c0)
                    +hy*(fp[-ml-1].c1-fp[-1].c1+fp[-ml].c1-fp->c1);
            }
        }
        // Solve the finite-element problem, and copy the pressure back into
        // the main data structure, subtracting off the mean, and taking into
        // account the boundary conditions
        mg2.solve_v_cycle();
        copy_pressure();

    // Update u and v based on c0, c1, and the computed pressure
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        while(fp<fe) {
            if(pres){
                fp->u=fp->c0-0.5*dt*rhoinv*xsp*(fp[ml+1].p+fp[1].p-fp[ml].p-fp->p);
                fp->v=fp->c1-0.5*dt*rhoinv*ysp*(fp[ml+1].p-fp[1].p+fp[ml].p-fp->p);
            }
            else{
                fp->u = fp->c0;
                fp->v = fp->c1;
            }
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
    double pavg=average_pressure();

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
    double pavg=0;

#pragma omp parallel for
    for(int j=0;j<n_fem;j++) {
        double *sop=sfem+j*m_fem,*soe=sop+(m_fem-1);
        double prow=*(sop++)*(x_prd?1:0.5);
        while(sop<soe) prow+=*(sop++);
        prow+=*sop*(x_prd?1:0.5);
        if(!y_prd&&(j==0||j==n)) prow*=0.5;
#pragma omp atomic
        pavg+=prow;
    }
    return pavg*(1./mn);
}

void fluid_2d::centered_diff(field *fp,double &uxx,double &uyy,double &vxx,double &vyy) {
    uyy=yysp*(fp[-ml].u-2*fp->u+fp[ml].u);
    vyy=yysp*(fp[-ml].v-2*fp->v+fp[ml].v);
    uxx=xxsp*(fp[-1].u-2*fp->u+fp[1].u);
    vxx=xxsp*(fp[-1].v-2*fp->v+fp[1].v);
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void fluid_2d::set_boundaries() {

    // Set left and right ghost values
    if(x_prd) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].u=fp[m-2].u;fp[-2].v=fp[m-2].v;
            fp[-1].u=fp[m-1].u;fp[-1].v=fp[m-1].v;
            fp[m].u=fp->u;fp[m].v=fp->v;
            fp[m+1].u=fp[1].u;fp[m+1].v=fp[1].v;
        }
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].u=-fp[1].u;fp[-2].v=-fp[1].v;
            fp[-1].u=-fp->u;fp[-1].v=-fp->v;
            fp[m].u=-fp[m-1].u;fp[m].v=-fp[m-1].v;
            fp[m+1].u=-fp[m-2].u;fp[m+1].v=-fp[m-2].v;
        }
    }

    // Set top and bottom ghost values
    const int tl=2*ml,g=n*ml;
    if(y_prd) {
        for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) {
            fp[-tl].u=fp[g-tl].u;fp[-tl].v=fp[g-tl].v;
            fp[-ml].u=fp[g-ml].u;fp[-ml].v=fp[g-ml].v;
            fp[g].u=fp->u;fp[g].v=fp->v;
            fp[g+ml].u=fp[ml].u;fp[g+ml].v=fp[ml].v;
        }
    } else {
        for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) {
            fp[-tl].u=-fp[ml].u;fp[-tl].v=-fp[ml].v;
            fp[-ml].u=-fp->u;fp[-ml].v=-fp->v;
            fp[g].u=-fp[g-ml].u;fp[g].v=-fp[g-ml].v;
            fp[g+ml].u=-fp[g-tl].u;fp[g+ml].v=-fp[g-tl].v;
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
            fp[-1].c0=fp[m-1].c0;fp[-1].c1=fp[m-1].c1;
        }
        xl=m;
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].c0=fp[-1].c1=0;
            fp[m].c0=fp[m].c1=0;
        }
        xl=m+1;
    }

    // Set top and bottom ghost values
    const int g=n*ml;
    if(y_prd) {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].c0=fp[g-ml].c0;fp[-ml].c1=fp[g-ml].c1;
        }
    } else {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].c0=fp[-ml].c1=0;
            fp[g].c0=fp[g].c1=0;
        }
    }
}

/** Sets up the fluid tracers by initializing them at random positions. */
void fluid_2d::init_tracers() {
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
    int i,j,i_prd,j_prd;
    double x,y,*tp=tm,*te=tm+(ntrace<<1);
    field *fp;
    while(tp<te) {

        // Find which grid cell the tracer is in
        //x=(*tp-ax)*xsp+0.5;y=(tp[1]-ay)*ysp+0.5; //potential bug
        x=(*tp-ax)*xsp-0.5;y=(tp[1]-ay)*ysp-0.5;
        i=floor(x);
        j=floor(y);

        // Compute the tracer's fractional position with the grid cell
        x-=i;y-=j;

        // Compute tracer's new position
        i_prd=x_prd?my_mod(i,m):i;
        j_prd=y_prd?my_mod(j,n):j;
        fp=fm+(i_prd+ml*j_prd);
        *(tp++)+=dt*((1-y)*(fp->u*(1-x)+fp[1].u*x)+y*(fp[ml].u*(1-x)+fp[ml+1].u*x));
        *(tp++)+=dt*((1-y)*(fp->v*(1-x)+fp[1].v*x)+y*(fp[ml].v*(1-x)+fp[ml+1].v*x));

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

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the tracer positions in batches of 128 floats
    int j,tbatch=(ntrace<<1)/128,tres=(ntrace<<1)%128;
    float *fp,*fe=buf+128;
    double *tp=tm,*te=tm+(ntrace<<1);

        double temp_x, temp_y;
    for(j=0;j<tbatch;j++) {
        fp=buf;
        while(fp<fe) {
                  temp_x = *(tp++);
                  temp_y = *(tp++);

                  remap_tracer(temp_x,temp_y);

                  *(fp++) = temp_x;
                  *(fp++) = temp_y;

                }
        fwrite(buf,sizeof(float),128,outf);
    }

    // Output the remaining tracer positions, if any
    if(tres>0) {
        fp=buf;
        do {
                  temp_x = *(tp++);
                  temp_y = *(tp++);

                  remap_tracer(temp_x,temp_y);

                  *(fp++) = temp_x;
                  *(fp++) = temp_y;

                } while(tp<te);
        fwrite(buf,sizeof(float),tres,outf);
    }

    // Close the file
    fclose(outf);
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void fluid_2d::write_files(const char *filename, int k) {
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/u.%d",filename, k);
    FILE *fu = safe_fopen(bufc, "w");
    sprintf(bufc,"%s/v.%d",filename, k);
    FILE *fv = safe_fopen(bufc, "w");
    sprintf(bufc, "%s/p.%d", filename, k);
    FILE *fp = safe_fopen(bufc, "w");
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            field *f = fm+j*ml+i;
            fprintf(fu, "%.8g ", f->u);
            fprintf(fv, "%.8g ", f->v);
            fprintf(fp, "%.8g ", f->p);
        }
        fprintf(fu, "\n");
        fprintf(fv, "\n");
        fprintf(fp, "\n");
    }
    fclose(fu);
    fclose(fv);
    fclose(fp);
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void fluid_2d::write_files(int k) {
    const int fflags=7;
    if(fflags&1) output("u",0,k);
    if(fflags&2) output("v",1,k);
    if(fflags&4) output("p",2,k);
    output_tracers("trace",k);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] ghost whether to output the ghost regions or not. */
void fluid_2d::output(const char *prefix,const int mode,const int sn,const bool ghost) {

    // Determine whether to output a cell-centered field or not
    bool cen=mode>=0&&mode<=1;
    int l=ghost?ml:(cen?m:m+1);
    double disp=(cen?0.5:0)-(ghost?2:0);

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    int i,j;
    float *bp=buf+1,*be=bp+l;
    *buf=l;
    for(i=0;i<l;i++) *(bp++)=ax+(i+disp)*dx;
    fwrite(buf,sizeof(float),l+1,outf);

    // Output the field values to the file
    //field *fr=ghost?fbase:fm;
    field *fr=ghost?fbase:err; //Plot error
    for(j=0;j<l;j++,fr+=ml) {
        field *fp=fr;
        *buf=ay+(j+disp)*dy;bp=buf+1;
        switch(mode) {
            case 0: while(bp<be) *(bp++)=(fp++)->u;break;
            case 1: while(bp<be) *(bp++)=(fp++)->v;break;
            case 2: while(bp<be) *(bp++)=(fp++)->p;break;
        }
        fwrite(buf,sizeof(float),l+1,outf);
    }

    // Close the file
    fclose(outf);
}
