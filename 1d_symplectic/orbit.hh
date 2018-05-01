#ifndef ORBIT_HH
#define ORBIT_HH

#include "sol_euler.hh"
#include "sol_improv_e.hh"
#include "sol_imp_onestep.hh"

#include <cstdio>
#include <cmath>

/** This class has functions to specify the test orbit problem. */
class orbit {
    public:
        double init_h;
        void orb_ff(double t_,double *in,double *out);
        void orb_init(double *q);
        void orb_print(double t_,double *q);
        double hamiltonian(double *q);
};

/** Class to solve the orbit problem with the forward Euler method. */
class orb_euler : public euler, public orbit {
    public:
        orb_euler() : euler(4) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_ff(t_,in,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the backward Euler method. */
class orb_back : public imp_onestep, public orbit {
    public:
        orb_back() : imp_onestep(4,1) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_ff(t_,in,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the improved Euler method. */
class orb_improv_e : public improv_e, public orbit {
    public:
        orb_improv_e() : improv_e(4) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_ff(t_,in,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the implicit midpoint method. */
class orb_imp_mid : public imp_onestep, public orbit {
    public:
        orb_imp_mid() : imp_onestep(4,0.5) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_ff(t_,in,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

#endif
