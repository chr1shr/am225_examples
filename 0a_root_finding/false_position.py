#!/usr/bin/python
from math import *
import sys

# Tolerance for convergence
tol=1e-14

# Function to consider
def f(x,lam):
    return lam-cos(x)

# Initial interval: assume f(a)<0 and f(b)>0
def false_position(lam):
    a=0
    b=pi
    fa=f(a,lam)
    fb=f(b,lam)

    # False position method
    it=0
    while it<100:

        # Evaluate function at the root of the linear interpolant
        c=a+(b-a)*fa/(fa-fb)
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
            a=c;fa=fc
        else:
            b=c;fb=fc
    print "# Too many iterations"
    sys.exit()

val=false_position(0.7)
print "# Root at",val
