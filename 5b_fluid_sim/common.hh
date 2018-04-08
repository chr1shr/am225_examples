#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cmath>

//#include <limits>
//const double big_number=std::numeric_limits<double>::max();

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);

#endif
