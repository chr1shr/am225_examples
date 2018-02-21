#include <cstdio>

int main() {

    // Variables must explicitly declared with a type
    int a=1,b;

    // Single-precision and double-precision floating point
    float c=2.0;
    double d=3.4;

    // Arithmetic
    b=3*(a++);

    // Formatted print
    printf("%d %d %g %g\n",a,b,c,d);
}
