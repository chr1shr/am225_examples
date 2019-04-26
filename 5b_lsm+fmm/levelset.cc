#include <cstring>
#include <limits>

#include "common.hh"
#include "levelset.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the level set.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] (x_prd_,y_prd_) the periodicity in the x and y directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] filename_ the filename of the output directory. */
levelset::levelset(const int m_,const int n_,const bool x_prd_,
        const bool y_prd_,const double ax_,const double bx_,
        const double ay_,const double by_,const char *filename_)
    : m(m_), n(n_), mn(m_*n_), ml(m+4), x_prd(x_prd_), y_prd(y_prd_), ax(ax_),
    ay(ay_), bx(bx_), by(by_), dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx),
    ysp(1/dy), xxsp(xsp*xsp), yysp(ysp*ysp), filename(filename_),
    fbase(new field[ml*(n+4)]), fm(fbase+2*ml+2),
    time(0.), f_num(0), oscillate_vel(false), buf(new float[m>123?m+5:128]) {}

/** The class destructor frees the dynamically allocated memory. */
levelset::~levelset() {
    delete [] buf;
    delete [] fbase;
}

/** Initializes the simulation, setting up the simulation fields, and choosing
 * the timestep.
 * \param[in] type the velocity field type.
 * \param[in] oscillate_vel_ whether to oscillate the velocity field.
 * \param[in] dt_pad_ the padding factor for the timestep, which should be
 *                    smaller than 1.
 * \param[in] max_spd a maximum speed from which to estimate the advection
 *                    timestep restriction. If a negative value is supplied,
 *                    then the advection CFL condition is explicitly
 *                    calculated. */
void levelset::initialize(int type,bool oscillate_vel_,double dt_pad,double max_spd) {
    oscillate_vel=oscillate_vel_;

    // Initialize the simulation fields
    init_fields(type);

    // Compute the timestep, based on the restrictions from the advection
    // and velocity, plus a padding factor
    choose_dt(dt_pad,max_spd<=0?advection_dt():(dx>dy?dx:dy)/max_spd);
}

/** Computes the maximum timestep that can resolve the advection, based on the
 * CFL condition.
 * \return The maximum timestep. */
double levelset::advection_dt() {
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

/** Chooses the timestep based on the limit from advection.
 * \param[in] dt_pad the padding factor for the timestep for the physical
 *                   terms, which should be smaller than 1.
 * \param[in] adv_dt the maximum timestep to resolve the advection.
 * \param[in] verbose whether to print out messages to the screen. */
void levelset::choose_dt(double dt_pad,double adv_dt,bool verbose) {
    dt_reg=dt_pad*adv_dt;

    // Print information if requested
    if(verbose) {
        printf("# Advection dt       : %g\n"
               "# Padding factor     : %g\n"
               "# Minimum dt         : %g\n",
               adv_dt,dt_pad,dt_reg);
    }
}

/** Initializes the simulation fields. */
void levelset::init_fields(int type) {

    // Loop over the primary grid and set the velocity and pressure
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay+dy*(j+0.5);
        field *fp=fm+ml*j;
        for(int i=0;i<m;i++,fp++) {
            double x=ax+dx*(i+0.5);
            switch(type) {
                case 0:
                    fp->u=y;fp->v=-x;break;
                case 1:
                    {
                        double r=sqrt(x*x+y*y),fr=sin(4*M_PI*r);
                        fp->u=fr*y;
                        fp->v=-fr*x;
                    }
            }
            double yy=y-0.5;
            if(yy<-1) yy+=2;
            fp->phi=sqrt(x*x+yy*yy)-0.3;
        }
    }

    // Now that the primary grid points are set up, initialize the ghost
    // points according to the boundary conditions
    set_boundaries();
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void levelset::solve(double duration,int frames) {
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
        printf("# Output frame %d [%d, %.8g s, %.8g s]\n",k,l,t1-t0,t2-t1);
        t0=t2;
    }
    f_num+=frames;
}

/** Steps the level set field forward.
 * \param[in] dt the time step to use. */
void levelset::step_forward(double dt) {
    int j;
    double hx=0.5*dt*xsp,hy=0.5*dt*ysp,
           ft=oscillate_vel?1:0.4*sin(time);

#pragma omp parallel for
    for(j=0;j<n;j++) for(int i=0;i<m;i++) {
        field *fp=fm+(ml*j+i),&f=*fp;

        // Compute advective terms using the second-order ENO scheme
        double phix,phiy,uc=ft*f.u,vc=ft*f.v;
        uc>0?vel_eno2(phix,hx,fp[1],f,fp[-1],fp[-2])
            :vel_eno2(phix,-hx,fp[-1],f,fp[1],fp[2]);
        vc>0?vel_eno2(phiy,hy,fp[ml],f,fp[-ml],fp[-2*ml])
            :vel_eno2(phiy,-hy,fp[-ml],f,fp[ml],fp[2*ml]);

        // Store the update to the level set field
        f.cphi=-uc*phix-vc*phiy;
    }

    // Apply the update
#pragma omp parallel for
    for(j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        while(fp<fe) (fp++)->update();
    }

    // Reset the ghost points according to the boundary conditions
    set_boundaries();
    time+=dt;
}

/** Calculates one-sided derivatives of the level set field using the
 * second-order ENO2 scheme.
 * \param[out] phid the computed ENO2 derivative.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with. */
inline void levelset::vel_eno2(double &phid,double hs,field &f0,field &f1,field &f2,field &f3) {
    phid=hs*eno2(f0.phi,f1.phi,f2.phi,f3.phi);
}

/** Calculates the ENO derivative using a sequence of values at four
 * gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double levelset::eno2(double p0,double p1,double p2,double p3) {
    return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void levelset::set_boundaries() {

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
            fp[-2].neu_bc(fp[1]);
            fp[-1].neu_bc(*fp);
            fp[m].neu_bc(fp[m-1]);
            fp[m+1].neu_bc(fp[m-2]);
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
            fp[-tl].neu_bc(fp[ml]);
            fp[-ml].neu_bc(*fp);
            fp[g].neu_bc(fp[g-ml]);
            fp[g+ml].neu_bc(fp[g-tl]);
        }
    }
}

/** Saves the header file.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void levelset::save_header(double duration, int frames) {
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/header",filename);
    FILE *outf=safe_fopen(bufc,f_num==0?"w":"a");
    fprintf(outf,"%g %g %d\n",time,time+duration,frames);
    fclose(outf);
}

/** Outputs the level set field to a file in a format that can be read by
 * Gnuplot.
 * \param[in] sn the current frame number to append to the filename.
 * \param[in] ghost whether to output the ghost regions or not. */
void levelset::output(const int sn,const bool ghost) {
    int l=ghost?ml:m;
    double disp=0.5-(ghost?2:0);

    // Assemble the output filename and open the output file
    char *bufc=((char*) buf);
    sprintf(bufc,"%s/phi.%d",filename,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    int i,j;
    float *bp=buf+1,*be=bp+l;
    *buf=l;
    for(i=0;i<l;i++) *(bp++)=ax+(i+disp)*dx;
    fwrite(buf,sizeof(float),l+1,outf);

    // Output the field values to the file
    field *fr=ghost?fbase:fm;
    for(j=0;j<l;j++,fr+=ml) {
        field *fp=fr;
        *buf=ay+(j+disp)*dy;bp=buf+1;
        while(bp<be) *(bp++)=(fp++)->phi;
        fwrite(buf,sizeof(float),l+1,outf);
    }

    // Close the file
    fclose(outf);
}
