#include <cstdio>

int main() {

    // Implement Fizz Buzz children's game
    for(int i=1;i<20;i++) {
        if(i%3==0) puts(i%5==0?"Fizz Buzz":"Fizz");
        else {
            if(i%5==0) puts("Buzz");
            else printf("%d\n",i);
        }
    }
}
