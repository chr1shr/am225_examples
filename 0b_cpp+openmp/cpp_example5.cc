#include <cstdio>

int main() {

    // Simple array construction
    double v[32];
    v[3]=4.;

    // A pointer to a double
    double* w;

    // Assign pointer. Note v itself is a pointer to the start
    // of the array.
    w=v+3;
    printf("%p %g\n",(void*) w,*w);

    // For-loop with pointers
    for(w=v;w<v+32;w++) *w=1.;

    // Out of bounds. May cause segmentation fault error. But may not.
    // With great power comes great responsibility.
    v[40]=0.;
}
