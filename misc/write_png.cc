#include "write_png.hh"

#include <cstdio>
#include <cstdlib>

#include <png.h>

inline png_byte wpng_truncate(double val,double amin,double amax) {
    return val<amin?0:(val>amax?255:static_cast<png_byte>(255.0/(amax-amin)*(val-amin)+0.5));
}

void wpng_abort(const char* err_msg) {
    fprintf(stderr,"PNG output: %s\n",err_msg);
    exit(1);
}

/** Writes a PNG image to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] (m,n) the dimensions of the image.
 * \param[in] co a 2D array of (R,G,B) color channels.
 * \param[in] (amin,amax) the values to scale to zero color and maximum color,
 *                        respectively. */
void write_png(const char *filename,int m,int n,double *co,double amin,double amax) {
    int j;

    // Convert
    png_bytep *rowp=new png_bytep[n];
    for(j=0;j<n;j++) {
        rowp[j]=new png_byte[3*m];
        for(png_byte* bp=rowp[j],*be=bp+3*m;bp<be;) {
            *(bp++)=wpng_truncate(*(co++),amin,amax);
            *(bp++)=wpng_truncate(*(co++),amin,amax);
            *(bp++)=wpng_truncate(*(co++),amin,amax);
        }
    }

    FILE *fp=fopen(filename,"wb");
    if(fp==NULL) wpng_abort("can't open output file");

    png_structp p=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    if(!p) wpng_abort("error creating write structure");

    png_infop info=png_create_info_struct(p);
    if(!p) wpng_abort("error creating info structure");

    if(setjmp(png_jmpbuf(p))) wpng_abort("setjmp error");
    png_init_io(p,fp);
    png_set_IHDR(p,info,m,n,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

    png_write_info(p,info);
    png_write_image(p,rowp);
    png_write_end(p,NULL);
    fclose(fp);

    // Free dynamically allocated memory and structures
    for(j=n-1;j>=0;j--) delete [] rowp[j];
    delete [] rowp;
    png_destroy_write_struct(&p,&info);
}
