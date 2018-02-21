#include <cstdio>

// Function to compute whether a number n is happy. Returns 0
// if happy, 1 if sad, and 2 if unresolved.
int happy(int n) {
    int iter=0,m,d;

    while(++iter<1000) {

        // Check for happiness (n=1) or sadness (n=4)
        if(n==1) return 0;
        if(n==4) return 1;

        // Map n to the next number in the sequence
        m=n;n=0;
        while(m>0) {

            // Extract lowest digit
            d=m%10;m/=10;

            // Add the square of the digit to the counter
            n+=d*d;
        }
    }

    // If the loop ended, then too many iterations were taken
    return 2;
}

int main() {

    // Specify dynamic threading: when a thread is done with its number, it
    // requests the next one available. Comes with performance overhead, but
    // good when the workload varies from number to number.
#pragma omp parallel for schedule(dynamic)
    for(int n=1;n<=100;n++) {
        switch(happy(n)) {
            case 0: printf("%d: that number is happy :)\n",n);break;
            case 1: printf("%d: that number is sad :(\n",n);break;
            case 2: printf("%d: that number is still in search of happiness :|\n",n);
        }
    }
}
