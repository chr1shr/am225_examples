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
    w(1), mem(2*(m+n)), he(new phi_field*[mem]) {
    setup_indicator_field();
}

/** The class destructor frees the dynamically allocated memory. */
fmm::~fmm() {
    delete [] buf;
    delete [] phibase;
}

/** Initializes the simulation fields. */
void fmm::init_fields() {
    const int jl=4*n/10,ju=6*n/10;

    for(int j=jl;j<ju;j++) {
        phi_field *fp=phim+ml*j+(m/2);

        fp->c=2;
        fp->phi=0;
    }
}

/** Initializes the heap. */
void fmm::init_heap() {
    w=0;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            phi_field *phip=phim+ml*j+i;
            if(phip->c==2) add_neighbors(phip);
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

double fmm::calc_phi(phi_field *phip) {
    double phiv,phih;
    bool vert=phi_look(phip,ml,phiv),
         horiz=phi_look(phip,1,phih);

    return vert?(horiz?phi_full(phiv,phih):phiv)
                  :(horiz?phih:0);
}

double fmm::phi_full(double phiv,double phih) {
    const double a=xxsp+yysp;
    double b=-phih*xxsp-phiv*yysp;
    double c=phih*phih*xxsp+phiv*phiv*yysp-1;
    return (-b+sqrt(b*b-a*c))/a;
}

inline bool fmm::phi_look(phi_field *phip,int d,double &phid) {
    if(phip[-d].c==2) {
        phid=phip[d].c==2?min(phip[-d].phi,phip[d].phi)
                         :phip[-d].c;
        return true;
    } else if(phip[d].c==2) {
        phid=phip[d].c;
        return true;
    }
    return false;
}

void fmm::add_to_heap(phi_field *phip) {
    double tphi=calc_phi(phip);
    int k=int(phip-phim),j=k/ml,i=k%ml;
    printf("%d %d %g\n",i,j,tphi);
	int c=w,bc=c>>1;
	while(bc>=1&&tphi<he[bc]->phi) {
		he[bc]->bp=c;
		he[c]=he[bc];
		c=bc;bc=c>>1;
	}
    he[c]=phip;
    he[c]->bp=c;
	he[c]->phi=tphi;
	he[c]->c=1;
	w++;
}

void fmm::trickle(phi_field *phip) {
    double tphi=calc_phi(phip);
	int c=phip->c,bc=c>>1;
	if(bc>=1&&tphi<he[bc]->phi) {
		do {
			he[bc]->bp=c;
			he[c]=he[bc];
			c=bc;bc=c>>1;
		} while(bc>=1&&tphi<he[bc]->phi);
        he[c]=phip;
        he[c]->bp=c;
	}
    he[c]->phi=tphi;
}

void fmm::fast_march() {
    phi_field *phip;
    init_heap();
    puts("yo");
    while(w>1) {
        phip=he[1];
        phip->c=2;
        reduce_heap();
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

void fmm::update(phi_field *phip) {
    if(phip->c==1) trickle(phip);
    else if(phip->c==0) add_to_heap(phip);
}

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

void fmm::add_heap_memory() {
    if(mem>m*n) {
        fputs("Maximum heap memory allocation exceeded\n",stderr);
        exit(1);
    }
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
