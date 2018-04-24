#include <cstring>
#include <limits>

#include "common.hh"
#include "fmm.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for level set.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds. */
fmm::fmm(const int m_,const int n_,const double ax_,const double bx_,
        const double ay_,const double by_)
    : m(m_), n(n_), mn(m_*n_), ml(m+4), ax(ax_), ay(ay_), bx(bx_), by(by_),
    dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx), ysp(1/dy), xxsp(xsp*xsp),
    yysp(ysp*ysp), phibase(new phi_field[ml*(n+4)]), phim(phibase+2*ml+2),
    w(1), mem(20*(m+n)), he(new phi_field*[mem]) {
    setup_indicator_field();
}

/** The class destructor frees the dynamically allocated memory. */
fmm::~fmm() {
    delete [] buf;
    delete [] phibase;
}

/** Initializes the simulation fields. */
void fmm::init_fields() {
    phi_field *fp=phim+10*ml+10;
    fp->c=2;
    fp->phi=0.;

/*    mark_box(45,55,30,80);
    return;

    const int jl=4*n/10,ju=6*n/10;
    for(int j=jl;j<ju;j++) {
        phi_field *fp=phim+ml*j+(m/2);
        fp->c=2;
        fp->phi=0;
    }*/
}

/** Marks a box of points to be part of the boundary. */
void fmm::mark_box(int il,int iu,int jl,int ju) {
    for(int j=jl;j<ju;j++) for(int i=il;i<iu;i++)
        phim[j*ml+i].c=3;
}

/** Initializes the heap by scanning over the grid and adding all the neighbors
 * of set points. */
void fmm::init_heap() {
    w=1;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            phi_field *phip=phim+ml*j+i;

            // If this point's indicator is 2, then scan all of its neighbors
            if(phip->c==2) add_neighbors(phip);
        }
    }
}

/** Scans the neighbors of a grid point, and adds any empty ones to the heap.
 * \param[in] phip a pointer to the grid point to consider. */
void fmm::add_neighbors(phi_field *phip) {

    // Check that there is enough heap memory for four new pointers
    if(w+4>mem) add_heap_memory();

    // Add any empty grid points to the heap
    if(phip[-ml].c==0) add_to_heap(phip-ml);
    if(phip[-1].c==0) add_to_heap(phip-1);
    if(phip[1].c==0) add_to_heap(phip+1);
    if(phip[ml].c==0) add_to_heap(phip+ml);
}

/** Calculates the value of the phi field at a given grid point.
 * \param[in] phip a pointer to the grid point to consider. */
double fmm::calc_phi(phi_field *phip) {
    double phiv,phih;

    // Look for available phi values in the horizontal and vertical directions,
    // finding the minimum
    bool vert=phi_look(phip,ml,phiv),
         horiz=phi_look(phip,1,phih);

    // Compute the phi value. For the case when both horizontal and vertical
    // neighbors are available, we need
    return vert?(horiz?phi_full(phiv,phih):phiv+dy)
               :(horiz?phih+dx:0);
}

double fmm::phi_full(double phiv,double phih) {
    const double a=xxsp+yysp;
    double b=-phih*xxsp-phiv*yysp;
    double c=phih*phih*xxsp+phiv*phiv*yysp-1;
    return (-b+sqrt(b*b-a*c))/a;
}

/** Finds the minimum set value of phi by comparing two adjacent neighbors.
 * \param[in] phip a pointer to the grid point to consider.
 * \param[in] d the memory step to the neighbors.
 * \param[out] phid the minimum set value, if available.
 * \return True if a value was found, false otherwise. */
inline bool fmm::phi_look(phi_field *phip,int d,double &phid) {
    if(phip[-d].c==2) {
        phid=phip[d].c==2?min(phip[-d].phi,phip[d].phi)
                         :phip[-d].phi;
        return true;
    } else if(phip[d].c==2) {
        phid=phip[d].phi;
        return true;
    }
    return false;
}

/** Adds a grid point to the heap.
 * \param[in] phip a pointer to the grid point to consider. */
void fmm::add_to_heap(phi_field *phip) {
    double tphi=calc_phi(phip);

    // Initially, the grid point will be placed at the top of the heap. Perform
    // swap operations until the heap property is restored.
    int c=w++,bc=c>>1;
    while(bc>=1&&tphi<he[bc]->phi) {
        he[bc]->bp=c;
        he[c]=he[bc];
        c=bc;bc=c>>1;
    }

    // Set the information for the new point
    he[c]=phip;
    he[c]->bp=c;
    he[c]->phi=tphi;
    he[c]->c=1;
}

/** Adjusts the heap when a phi value may have decreased. */
void fmm::trickle(phi_field *phip) {

    // Compute the new phi value
    double tphi=calc_phi(phip);

    // If the phi value is smaller than the parent, then swap it with the parent
    int c=phip->bp,bc=c>>1;
    if(bc>=1&&tphi<he[bc]->phi) {
        do {
            he[bc]->bp=c;
            he[c]=he[bc];
            c=bc;bc=c>>1;
        } while(bc>=1&&tphi<he[bc]->phi);
        he[c]=phip;
        he[c]->bp=c;
    }

    // Update the new phi value
    he[c]->phi=tphi;
}

/** Performs the fast marching method. */
void fmm::fast_march() {
    phi_field *phip;

    // Initialize the heap, by finding neighbors of set points
    init_heap();

    // Loop until the heap is non-empty, each time considering the point with
    // the smallest phi value
    while(w>1) {

        // Change the status of the point with the smallest phi value to "set"
        phip=he[1];
        phip->c=2;

        // Remove this point from the heap
        reduce_heap();

        // Check neighbors of this point, to see if they must be added to the
        // heap or have their phi values updated
        if(w+4>mem) add_heap_memory();
        update(phip-ml);update(phip-1);
        update(phip+1);update(phip+ml);
    }
}

void fmm::reduce_heap() {
    phi_field *phip=he[--w];
    double &tphi=phip->phi;
    int bc=1,c=bc<<1,cmin;
    while(c+1<w) {
        if(he[c]->phi<tphi) {
            cmin=he[c+1]->phi<he[c]->phi?c+1:c;
        } else {
            if(he[c+1]->phi<tphi) cmin=c+1;
            else break;
        }
        he[bc]=he[cmin];
        he[bc]->bp=bc;
        bc=cmin;
        c=bc<<1;
    }
    if(c+1==w) {
        if(he[c]->phi<tphi) {
            he[bc]=he[c];
            he[bc]->bp=bc;
            bc=c;
        }
    }
    he[bc]=phip;
    he[bc]->bp=bc;
}

/* Updates a grid point during the fast marching calculation. If it has never been
 * considered, it is added to the heap. If it is already marked, then its phi value
 * is updated.
 * \param[in] phip a pointer to the grid point to consider. */
void fmm::update(phi_field *phip) {
    if(phip->c==1) trickle(phip);
    else if(phip->c==0) add_to_heap(phip);
}

/** Sets up the indicator field that marks the status of the grid points during
 * the fast march. */
void fmm::setup_indicator_field() {
    phi_field *phip=phibase,*phie;
    for(phie=phip+2*ml;phip<phie;phip++) phip->c=3;
    while(phip<phibase+(n+2)*ml) {
        (phip++)->c=3;(phip++)->c=3;
        for(phie=phip+m;phip<phie;phip++) phip->c=0;
        (phip++)->c=3;(phip++)->c=3;
    }
    for(phie=phip+2*ml;phip<phie;phip++) phip->c=3;
}

/** Doubles the size of the heap memory. */
void fmm::add_heap_memory() {

    // If the current allocation exceeds the total number of grid points, that indicates
    // a problem, so exit with an error.
    if(mem>m*n) {
        fputs("Maximum heap memory allocation exceeded\n",stderr);
        exit(1);
    }

    // Allocate an array of double the size, and copy the contents of the
    // current array into it
    mem<<=1;
    phi_field **nhe=new phi_field*[mem];
    for(int i=0;i<w;i++) nhe[i]=he[i];
    delete [] he;
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] filename the field name to use as the filename prefix. */
void fmm::output(const char *filename,int mode) {

    // Assemble the output filename and open the output file
    FILE *outf=safe_fopen(filename,"wb");

    // Output the first line of the file
    int i,j;
    float *buf=new float[m+1],*bp=buf+1,*be=bp+m;
    *buf=m;
    for(i=0;i<m;i++) *(bp++)=ax+(i+0.5)*dx;
    fwrite(buf,sizeof(float),m+1,outf);

    // Output the field values to the file
    phi_field *phir=phim;
    for(j=0;j<m;j++,phir+=ml) {
        phi_field *phip=phir;
        *buf=ay+(j+0.5)*dy;bp=buf+1;
        switch(mode) {
            case 0: while(bp<be) *(bp++)=(phip++)->phi;break;
            case 1: while(bp<be) *(bp++)=(phip++)->c;
        }
        fwrite(buf,sizeof(float),m+1,outf);
    }

    // Close the file and free the dynamically allocated memory
    fclose(outf);
    delete [] buf;
}

/** Sets the boundary points to have a large phi value. */
void fmm::set_boundary_phi() {
    const double phimax=1.05*sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
    for(int i=0;i<ml*(n+4);i++) if(phibase[i].c==3) phibase[i].phi=phimax;
}

/** Integrates an ODE following the negative gradient.
 * \param[in] (x,y) the starting point. */
void fmm::integrate_path(double x,double y) {
    double dt=0.0001;
    double p,gx,gy;

    p=bilinear(x,y,gx,gy);
    printf("%g %g %g %g %g\n",x,y,p,gx,gy);

    while(p>0.1*(dx+dy)) {
        x-=dt*gx;
        y-=dt*gy;
        p=bilinear(x,y,gx,gy);
        printf("%g %g %g %g %g\n",x,y,p,gx,gy);
    }
}

double fmm::bilinear(double x,double y,double &gx,double &gy) {
    x=(x-ax)*xsp-0.5;y=(y-ay)*ysp-0.5;
    int i=int(x);
    int j=int(y);
    if(i<-1) i=-1;else if(i>=m) i=m-1;
    if(j<-1) j=-1;else if(j>=n) i=n-1;

    // Compute the tracer's fractional position with the grid cell
    x-=i;y-=j;

    // Compute tracer's new position
    phi_field *fp=phim+(i+ml*j);
    gx=xsp*((1-y)*(fp[1].phi-fp->phi)+y*(fp[ml+1].phi-fp[ml].phi));
    gy=ysp*((1-x)*(fp[ml].phi-fp->phi)+x*(fp[ml+1].phi-fp[1].phi));
    return (1-y)*(fp->phi*(1-x)+fp[1].phi*x)+y*(fp[ml].phi*(1-x)+fp[ml+1].phi*x);
}
