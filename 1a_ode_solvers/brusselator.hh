#ifndef BRUSSELATOR_HH
#define BRUSSELATOR_HH

#include "sol_euler.hh"
#include "sol_ralston.hh"
#include "sol_heun3.hh"
#include "sol_rk4.hh"
#include "sol_hammer_h.hh"

/** This class has functions to specify the test Brusselator problem. */
class brusselator {
    public:
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void brus_ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void brus_init(double *q) {
            *q=1.5;
            q[1]=3.;
        }
};

/** Class to solve the Brusselator problem with the Euler method. */
class brus_euler : public euler, public brusselator {
    public:
        brus_euler() : euler(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

/** Class to solve the Brusselator problem with the Ralston method. */
class brus_ralston : public ralston, public brusselator {
    public:
        brus_ralston() : ralston(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

/** Class to solve the Brusselator problem with the third-order Heun method. */
class brus_heun3 : public heun3, public brusselator {
    public:
        brus_heun3() : heun3(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

/** Class to solve the Brusselator problem with the fourth-order Runge-Kutta
 * method. */
class brus_rk4 : public rk4, public brusselator {
    public:
        brus_rk4() : rk4(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

/** Class to solve the Brusselator problem with the fourth-order
 * Hammer-Hollingsworth method. */
class brus_hammer_h : public hammer_h, public brusselator {
    public:
        brus_hammer_h() : hammer_h(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

#endif
