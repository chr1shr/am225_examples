#include <cstring>
#include <limits>

#include "common.hh"
#include "fmm.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for level set.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] (x_prd_,y_prd_) the periodicity in the x and y directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds. */
fmm::fmm(const int m_,const int n_,const bool x_prd_,
        const bool y_prd_,const double ax_,const double bx_,
        const double ay_,const double by_,const char *filename_)
    : m(m_), n(n_), mn(m_*n_), ml(m+4), x_prd(x_prd_), y_prd(y_prd_), ax(ax_),
    ay(ay_), bx(bx_), by(by_), dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx),
    ysp(1/dy), phibase(new phi_field[ml*(n+4)]),
    phim(fbase+2*ml+2), w(0), mem(2*(m+n)), he(new phi_field*[mem]) {

    // Initialize the indicator field, and set the ghost regions
    for(int i=0;i<ml*(n+4)) phibase[i].c=0;

    for(int i=0;i<ml;
}

/** The class destructor frees the dynamically allocated memory. */
fmm::~fmm() {
    delete [] buf;
    delete [] phibase;
}

/** Initializes the simulation fields. */
void fmm::init_fields(int type) {
    const int jl=4*n/10,ju=6*n/10;

    for(int j=jl;j<ju;j++) {
        phi_field *fp=phim+ml*j+(m/2);

        fp->c=2;
        fp->phi=0;
    }
}

void fmm::init_heap() {
    w=0;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            phi_field *phip=phim+ml*j+i;
            if(*phip==2) add_neighbors(*phip);
        }
    }
}

void fmm::add_neighbors(phi_field *phip) {
    if(w+4>mem) add_heap_memory();
    if(phip[-ml].c==0) add_to_heap(phip-ml);
    if(phip[-1].c==0) add_to_heap(phip-1);
    if(phip[1].c==0) add_to_heap(phip+1);
    if(phip[ml].c==0) add_to_heap(phip+ml);
}

void fmm::add_to_heap(phi_field *phip) {
    *(hep++)=phip;
    calc_phi(phip);
}

void fmm:calc_phi(phi_field *phip) {
    if(phip[-1].c&3==2) {
}

void fmm::introduce(int ij,double tphi) {
	int c=w,bc=c>>1;
	while(bc>=1&&tphi<he[bc].phi) {
		he[bc].bp=c;
		he[c]=he[bc];
		c=bc;bc=c>>1;
	}
    he[c]=phim[ij];
    he[c].bp=c;
	he[c].phi=tphi;
	hp[c].c=1;
	w++;
}

void fmm::trickle(int c,int ij,double tphi) {
	int bc=c>>1;
	if(bc>=1&&tphi<he[bc].phi) {
		do {
			he[bc].bp=c;
			he[c]=he[bc];
			c=bc;bc=c>>1;
		} while(bc>=1&&tphi<he[bc].phi);
		
        he[c]=bp.c;
        bp[ij]=c;
		hp[c]=ij;
	}
	phi[ij]=tphi;
}

void fmm::update(int i,int ij) {
	if (s[ij]==z.q) {
		s[ij]=z.p;
		if(w==mem) add_memory();
		introduce(ij,calc(i,ij));
	} else if (s[ij]==z.p) {
		trickle(bp[ij],ij,calc(i,ij));
	}
}

void fmm::set(int i,int ij) {

}


void fmm::add_heap_memory() {
    mem>>=1;
    if(mem>=max_heap_memory) {
        fputs("Maximum heap memory allocation exceeded\n",stderr);
        exit(1);
    }
    nhe=new phi_field*[mem];
    for(int i=0;i<w;i++) nhe[i]=he[i];
    delete [] he;
}

