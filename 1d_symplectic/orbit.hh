#ifndef ORBIT_HH
#define ORBIT_HH

#include "sol_euler.hh"
#include "sol_improv_e.hh"
#include "sol_imp_onestep.hh"
#include "sol_sep.hh"

#include <cstdio>
#include <cmath>

/** This class has functions to specify the test orbit problem. */
class orbit {
    public:
        /** The semi-major axis of the orbit. */
        const double a;
        /** The eccentricity of the orbit. */
        const double e;
        orbit(double a_,double e_) : a(a_), e(e_) {}
        /** The initial value of the Hamiltonian, which is used to track
         * changes in the value over time. */
        double init_h;
        /** Computes the change in the p (momentum) variables from the q
         * (coordinate) variables.
         * \param[in] in the current state of the q variables.
         * \param[in] out the change in the p variables. */
        inline void orb_fq(double *in,double *out) {
            *out=*in;
            out[1]=in[1];
        }
        /** Computes the change in the q (momentum) variables from the p
         * (coordinate) variables.
         * \param[in] in the current state of the p variables.
         * \param[in] out the change in the q variables. */
        inline void orb_fp(double *in,double *out) {
            double fac=1/(*in*(*in)+in[1]*in[1]);
            fac*=sqrt(fac);
            *out=-*in*fac;
            out[1]=-in[1]*fac;
        }
        void orb_init(double *q);
        void orb_print(double t_,double *q);
        double hamiltonian(double *q);
};

/** Class to solve the orbit problem with the forward Euler method. */
class orb_euler : public euler, public orbit {
    public:
        orb_euler(double a_,double e_) : euler(4), orbit(a_,e_) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_fp(in,out+2);orb_fq(in+2,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the backward Euler method. */
class orb_back : public imp_onestep, public orbit {
    public:
        orb_back(double a_,double e_) : imp_onestep(4,1), orbit(a_,e_) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_fp(in,out+2);orb_fq(in+2,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the improved Euler method. */
class orb_improv_e : public improv_e, public orbit {
    public:
        orb_improv_e(double a_,double e_) : improv_e(4), orbit(a_,e_) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_fp(in,out+2);orb_fq(in+2,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the implicit midpoint method. */
class orb_imp_mid : public imp_onestep, public orbit {
    public:
        orb_imp_mid(double a_,double e_) : imp_onestep(4,0.5), orbit(a_,e_) {}
        virtual void ff(double t_,double *in,double *out) {
            orb_fp(in,out+2);orb_fq(in+2,out);
        }
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with the first-order symplectic method. */
class orb_fo_symplectic : public sym_separable, public orbit {
    public:
        orb_fo_symplectic(double a_,double e_) : sym_separable(4,0), orbit(a_,e_) {}
        virtual void fp(double *in,double *out) {orb_fp(in,out);}
        virtual void fq(double *in,double *out) {orb_fq(in,out);}
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with Ruth's third-order method. */
class orb_ruth3 : public sym_separable, public orbit {
    public:
        orb_ruth3(double a_,double e_) : sym_separable(4,1), orbit(a_,e_) {}
        virtual void fp(double *in,double *out) {orb_fp(in,out);}
        virtual void fq(double *in,double *out) {orb_fq(in,out);}
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

/** Class to solve the orbit problem with Ruth's fourth-order method (based on
 * concatenation with the adjoint). */
class orb_ruth4 : public sym_separable, public orbit {
    public:
        orb_ruth4(double a_,double e_) : sym_separable(4,2), orbit(a_,e_) {}
        virtual void fp(double *in,double *out) {orb_fp(in,out);}
        virtual void fq(double *in,double *out) {orb_fq(in,out);}
        virtual void print() {orb_print(t,q);}
        virtual void init() {orb_init(q);}
};

#endif
