#!/usr/bin/python
from math import *
import sys

# Tolerance for convergence
tol=1e-14

# Function to consider
def f(x,lam):
    return lam-cos(x)

# Initial interval: assume f(a)<0 and f(b)>0
def bisection(lam):
    a=0
    b=pi

    # Bisection search
    it=0
    while it<100:

        # Midpoint evaluation
        c=0.5*(b+a)
        fc=f(c,lam)

        # Check for acceptable solution
        if it>0 and abs(prev_c-c)<tol:
            return c

        # Print status message
        print it,a,b,c,c-acos(lam)
        it+=1
        prev_c=c

        # Reduce interval size
        if fc<0:
            a=c
        else:
            b=c
    print "# Too many iterations"
    sys.exit()

val=bisection(0.7)
print "# Root at",val
