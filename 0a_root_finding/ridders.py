#!/usr/bin/python
from math import *
import sys

# Tolerance for convergence
tol=1e-14

# Function to consider
def f(x,lam):
    return lam-cos(x)

# Initial interval: assume f(a)<0 and f(b)>0
def ridders(lam):
    a=0
    b=pi
    fa=f(a,lam)
    fb=f(b,lam)

    # Ridders root finding
    it=0
    while it<100:

        # Midpoint evaluation
        c=0.5*(b+a)
        fc=f(c,lam)

        # Ridders calculation
        d=c-(c-a)*((fc,-fc)[fa>fb])/sqrt(fc*fc-fa*fb)
        fd=f(d,lam)

        # Check for acceptable solution
        if it>0 and abs(prev_d-d)<tol:
            return d

        # Print status message
        print it,a,b,d,d-acos(lam)
        it+=1
        prev_d=d

        # Reduce interval size
        if fd<0:
            a=d;fa=fd
            if fc>=0:
                b=c;fb=fc
        else:
            b=d;fb=fd
            if fc<0:
                a=c;fa=fc
    print "# Too many iterations"
    sys.exit()

val=ridders(0.7)
print "# Root at",val
