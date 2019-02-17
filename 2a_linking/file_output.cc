#include "file_output.hh"

/** Outputs a two dimensional array in the Gnuplot binary format.
 * \param[in] filename the name of the file to write to.
 * \param[in] fld the array to use.
 * \param[in] m the horizontal grid size.
 * \param[in] n the vertical grid size. */
void gnuplot_output(const char* filename,double *fld,int m,int n,double ax,double bx,double ay,double by) {
    double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);
    int i,j;

    // Open file and print error if there is a problem
    FILE *outf=fopen(filename,"wb");
    if(outf==NULL) {
        fputs("Can't open file\n",stderr);
        exit(1);
    }

    // Allocate memory and write the header file
    float *fbuf=new float[m+1],*fp=fbuf;
    double *pp;
    *(fp++)=m;for(i=0;i<m;i++) *(fp++)=ax+i*dx;
    sfwrite(fbuf,sizeof(float),m+1,outf);
    
    // Write field entries line-by-line
    for(j=0;j<n;j++) {

        // Write header entry
        fp=fbuf;*(fp++)=ay+j*dy;

        // Write a horizontal line to the buffer
        pp=fld+j*m;for(i=0;i<m;i++) *(fp++)=static_cast<float>(*(pp++));
        sfwrite(fbuf,sizeof(float),m+1,outf);
    }

    // Remove temporary memory and close file
    delete [] fbuf;
    fclose(outf);
}
