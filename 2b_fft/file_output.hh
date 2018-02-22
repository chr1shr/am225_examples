#ifndef FILE_OUTPUT_HH
#define FILE_OUTPUT_HH

#include <cstdio>
#include <cstdlib>

/** Performs an fwrite operation, and prints out an error in the case
 * when it fails.
 * \param[in] bu a pointer to the location to write.
 * \param[in] size the size in bytes of the entries to write.
 * \param[in] count the number of entries to write.
 * \param[in] f the file handle to write to. */
inline void sfwrite(const void *bu,size_t size,size_t count,FILE *f) {
    if(fwrite(bu,size,count,f)!=count) {
        fprintf(stderr,"Error writing data to file");
        exit(1);
    }
}

void gnuplot_output(const char* filename,double *fld,int m,int n,double ax,double bx,double ay,double by);

#endif
