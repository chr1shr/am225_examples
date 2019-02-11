#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <cstdio>
#include <cstdlib>

#include "dop853.hh"

/** Opens a file, and checks the return value to ensure that the operation was
 * successful.
 * \param[in] filename the file to open.
 * \return The file handle. */
FILE* safe_fopen(const char *filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fprintf(stderr,"Can't open output file '%s'\n",filename);
        exit(1);
    }
    return fp;
}

/** The Brusselator test problem from Hairer et al. */
class brusselator : public dop853 {
    public:
        FILE *fp;
        brusselator() : dop853(2) {
            fp=safe_fopen("dt.out1");
                *w=1.5;w[1]=3.;
        }
        ~brusselator() {
            fclose(fp);
        }
        virtual void ff(double tt,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        virtual void output() {
            fprintf(fp,"%.15g %.15g %.15g %.15g\n",
                t,*w,w[1],h);
        }
};

#endif
