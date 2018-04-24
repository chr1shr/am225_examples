#ifndef WRITE_PNG_HH
#define WRITE_PNG_HH

void wpng_abort(const char* err_msg);
void write_png(const char *filename,int m,int n,double *co,double amin,double amax);

#endif
