// Advection-diffusion-limited dissolution code
//
// Author   : Chris H. Rycroft (Harvard SEAS / Lawrence Berkeley Laboratory)
// Email    : chr@alum.mit.edu
// Date     : May 6th 2015
//
// This source code file is based on the adaptive eighth-order Dormand-Prince
// integration code DOP853.F, developed by E. Hairer and G. Wanner. The
// original code is written is Fortran but has been refactored here into a C++
// class. For more information, see
//
// http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
//
// E. Hairer, S.P. Norsett AND G. Wanner, Solving ordinary differential
// equations I. Nonstiff problems, 2nd edition. Springer series in
// computational mathematics, Springer-Verlag (1993).

/** \file dop853.cc
 * \brief Function implementations for the dop853 class. */

#include <cstdio>
#include <cmath>
#include <cstring>

#include "dop853.hh"

/** The construct allocates the memory for the intermediate steps during the
 * time-integration, and sets constants that are used in the method.
 * \param[in] n_ the number of degrees of freedom in the problem. */
dop853::dop853(int n_) : n(n_), nmax(16777216), nstiff(4), fac1(1/3.0),
	fac2(6.0), safe(0.9), beta(0.0), uround(2.3e-16), atoli(1e-14),
	rtoli(0), w(new double[n]), ww1(new double[n]), k1(new double[n]),
	k2(new double[n]), k3(new double[n]), k4(new double[n]),
	k5(new double[n]), k6(new double[n]), k7(new double[n]),
	k8(new double[n]), k9(new double[n]), k10(new double[n]),
	rc1(new double[n]), rc2(new double[n]), rc3(new double[n]),
	rc4(new double[n]), rc5(new double[n]), rc6(new double[n]),
	rc7(new double[n]), rc8(new double[n]) {}

/** The destructor frees the dynamically allocated memory. */
dop853::~dop853() {
	delete [] rc8;delete [] rc7;delete [] rc6;delete [] rc5;
	delete [] rc4;delete [] rc3;delete [] rc2;delete [] rc1;
	delete [] k10;delete [] k9;delete [] k8;delete [] k7;
	delete [] k6;delete [] k5;delete [] k4;delete [] k3;
	delete [] k2;delete [] k1;delete [] ww1;delete [] w;
}

/** Sets the counters that are used for diagnostic purposes. */
void dop853::init_counters() {
	nstep=0;
	naccpt=0;
	nrejct=0;
	nfixed=0;
	nff=0;
	nonsti=0;
	iasti=0;
}

/** Solves the ODE problem over a given time interval, storing snapshots
 * of the fields at specific times.
 * \param[in] (tstart,tend) the time interval to integrate over.
 * \param[in] snaps the number of snapshots to store.
 * \param[in] ws an array of pointers to where to store the snapshots. */
int dop853::solve(double tstart,double tend,int snaps,double **ws) {
	const double facc1=1/fac1,facc2=1/fac2,expo1=1.0/8-beta*0.2;
	double fac11,fac,posneg=tend>tstart?1:-1,facold=1e-4,hmax=fabs(tend-tstart),
	       hnew,err,sf=snaps>1?(tend-tstart)/(snaps-1):0,isf=snaps>1?1.0/sf:0;
	bool last=false,reject=false;
	t=0;csn=0;int nsn;
	init_counters();

	// Do any problem-specific output and save the first snapshot if needed
	output();
	if(snaps>0) memcpy(*ws,w,n*sizeof(double));

	// Do initial function call, and use it to estimate the initial step
	// size
	ff(t,w,k1);nff++;
	h=hinit(hmax,posneg);
	while(true) {

		// Check for the case when too many integration steps have been
		// perform
		if(nstep>nmax) {
			if(snaps>1) memcpy(ws[++csn],w,n*sizeof(double));
			return 2;
		}

		// Check for the case when the timestep has become too small
		if(0.1*fabs(h)<=fabs(t)*uround) {
			if(snaps>1) memcpy(ws[++csn],w,n*sizeof(double));
			return 3;
		}

		// Check for the case when the integration is close to
		// completion, and if so, try stepping to the end of the
		// interval
		if((t+1.01*h-tend)*posneg>0.0) {
			h=tend-t;
			last=true;
    		}

		// Perform the eighth-order integration step, and estimate the
		// local error
		nstep++;
		step12(k5);
		err=fabs(h)*error_estimation();

		// Estimate new timestep based on the local error
		fac11=pow(err,expo1);
		fac=fac11*std::pow(facold,-beta);
		fac=max_d(facc2,min_d(facc1,fac/safe));
		hnew=h/fac;

		// Check whether the estimated error is within the tolerance
		if(err<=1.0) {

			// If the error is within the tolerance, then accept
			// the step
			facold=max_d(err,1.0E-4);
			naccpt++;
			ff(tph,k5,k4);
			nff++;

			// Carry out stiffness detection at periodic intervals,
			// or if stiffness has previously been detected
			if((!(naccpt%nstiff)||(iasti>0))&&detect_stiffness()) {
				if(snaps>1) memcpy(ws[++csn],w,n*sizeof(double));
				return 1;
			}

			// Check whether snapshots are required over this integration
			// step, and if so, compute the dense output formulae
			nsn=int((tph-tstart)*isf);
			if(nsn>=snaps-1) nsn=snaps-2;
			if(nsn>csn) dense_output();

			// Copy the computed step into the state vector, and
			// the last Runge--Kutta step into the new first
			// Runge--Kutta step. Update the time.
			memcpy(k1,k4,n*sizeof(double));
			memcpy(w,k5,n*sizeof(double));
			told=t;
			t=tph;

			// Carry out any timestep-specific output
			output();

			// Compute any snapshots using the dense output
			while(csn<nsn) {csn++;dense(ws[csn],tstart+csn*sf);}

			// If this is the last timestep, then store a snapshot
			// and return
			if(last) {
				if(snaps>1) memcpy(ws[++csn],w,n*sizeof(double));
				return 0;
			}

			// Check for special cases for timestep choice
			if(fabs(hnew)>hmax) hnew=posneg*hmax;
			if(reject) hnew=posneg*min_d(fabs(hnew),fabs(h));
			reject=false;
		} else {

			// If the local error exceeded the tolerance, then
			// compute a new timestep and try again. Note that as
			// in the original DOP853.F, rejected steps at the
			// start of the computation are not counted.
			hnew=h/min_d(facc1,fac11/safe);
			reject=true;
			if(naccpt>=1) nrejct++;
			last=false;
		}
		h=hnew;
	}
}

/** Takes a step of a fixed size.
 * \param[in] h_ the step size to take.
 * \param[in] last whether this is the last timestep or not. */
void dop853::fixed_step(double h_,bool last) {
	h=h_;
	step12(w);
	t=tph;
	if(!last) {ff(t,w,k1);nff++;}
	nstep++;nfixed++;
}

/** Carries out the twelve steps of the eighth-order Dormand-Prince scheme.
 * \param[in] p1 a pointer to where to store the final calculation step. */
void dop853::step12(double *p1) {
	int i;
	const double    c2  = 0.526001519587677318785587544488E-01,
			c3  = 0.789002279381515978178381316732E-01,
			c4  = 0.118350341907227396726757197510E+00,
			c5  = 0.281649658092772603273242802490E+00,
			c6  = 0.333333333333333333333333333333E+00,
			c7  = 0.25E+00,
			c8  = 0.307692307692307692307692307692E+00,
			c9  = 0.651282051282051282051282051282E+00,
			c10 = 0.6E+00,
			c11 = 0.857142857142857142857142857142E+00,

			b1 =   5.42937341165687622380535766363E-2,
			b6 =   4.45031289275240888144113950566E0,
			b7 =   1.89151789931450038304281599044E0,
			b8 =  -5.8012039600105847814672114227E0,
			b9 =   3.1116436695781989440891606237E-1,
			b10 = -1.52160949662516078556178806805E-1,
			b11 =  2.01365400804030348374776537501E-1,
			b12 =  4.47106157277725905176885569043E-2,
			a21 =    5.26001519587677318785587544488E-2,
			a31 =    1.97250569845378994544595329183E-2,
			a32 =    5.91751709536136983633785987549E-2,
			a41 =    2.95875854768068491816892993775E-2,
			a43 =    8.87627564304205475450678981324E-2,
			a51 =    2.41365134159266685502369798665E-1,
			a53 =   -8.84549479328286085344864962717E-1,
			a54 =    9.24834003261792003115737966543E-1,
			a61 =    3.7037037037037037037037037037E-2,
			a64 =    1.70828608729473871279604482173E-1,
			a65 =    1.25467687566822425016691814123E-1,
			a71 =    3.7109375E-2,
			a74 =    1.70252211019544039314978060272E-1,
			a75 =    6.02165389804559606850219397283E-2,
			a76 =   -1.7578125E-2,

			a81 =    3.70920001185047927108779319836E-2,
			a84 =    1.70383925712239993810214054705E-1,
			a85 =    1.07262030446373284651809199168E-1,
			a86 =   -1.53194377486244017527936158236E-2,
			a87 =    8.27378916381402288758473766002E-3,
			a91 =    6.24110958716075717114429577812E-1,
			a94 =   -3.36089262944694129406857109825E0,
			a95 =   -8.68219346841726006818189891453E-1,
			a96 =    2.75920996994467083049415600797E1,
			a97 =    2.01540675504778934086186788979E1,
			a98 =   -4.34898841810699588477366255144E1,
			a101 =   4.77662536438264365890433908527E-1,
			a104 =  -2.48811461997166764192642586468E0,
			a105 =  -5.90290826836842996371446475743E-1,
			a106 =   2.12300514481811942347288949897E1,
			a107 =   1.52792336328824235832596922938E1,
			a108 =  -3.32882109689848629194453265587E1,
			a109 =  -2.03312017085086261358222928593E-2,

			a111 =  -9.3714243008598732571704021658E-1,
			a114 =   5.18637242884406370830023853209E0,
			a115 =   1.09143734899672957818500254654E0,
			a116 =  -8.14978701074692612513997267357E0,
			a117 =  -1.85200656599969598641566180701E1,
			a118 =   2.27394870993505042818970056734E1,
			a119 =   2.49360555267965238987089396762E0,
			a1110 = -3.0467644718982195003823669022E0,
			a121 =   2.27331014751653820792359768449E0,
			a124 =  -1.05344954667372501984066689879E1,
			a125 =  -2.00087205822486249909675718444E0,
			a126 =  -1.79589318631187989172765950534E1,
			a127 =   2.79488845294199600508499808837E1,
			a128 =  -2.85899827713502369474065508674E0,
			a129 =  -8.87285693353062954433549289258E0,
			a1210 =  1.23605671757943030647266201528E1,
			a1211 =  6.43392746015763530355970484046E-1;

	for(i=0;i<n;i++) ww1[i]=w[i]+h*a21*k1[i];
	ff(t+c2*h,ww1,k2);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a31*k1[i]+a32*k2[i]);
	ff(t+c3*h,ww1,k3);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a41*k1[i]+a43*k3[i]);
	ff(t+c4*h,ww1,k4);
	for(i=0;i <n;i++) ww1[i]=w[i]+h*(a51*k1[i]+a53*k3[i]+a54*k4[i]);
	ff(t+c5*h,ww1,k5);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a61*k1[i]+a64*k4[i]+a65*k5[i]);
	ff(t+c6*h,ww1,k6);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a71*k1[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
	ff(t+c7*h,ww1,k7);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a81*k1[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
	ff(t+c8*h,ww1,k8);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);
	ff(t+c9*h,ww1,k9);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a101*k1[i]+a104*k4[i]+a105*k5[i]
			       +a106*k6[i]+a107*k7[i]+a108*k8[i]+a109*k9[i]);
	ff(t+c10*h,ww1,k10);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a111*k1[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]
			       +a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
	ff(t+c11*h,ww1,k2);
	tph=t+h;
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a121*k1[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]
			       +a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k2[i]);
	ff(tph,ww1,k3);
	nff+=11;
	for(i=0;i<n;i++) {
		k4[i]=b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k2[i]+b12*k3[i];
		p1[i]=w[i]+h*k4[i];
	}
}

/** Estimates the local numerical error, using a combination of a fifth-order
 * formula and a third-order formula.
 * \return The calculated error. */
double dop853::error_estimation() {
	const double bhh1=0.244094488188976377952755905512E+00,
		     bhh2=0.733846688281611857341361741547E+00,
		     bhh3=0.220588235294117647058823529412E-01,
		     er1 = 0.1312004499419488073250102996E-01,
		     er6 =-0.1225156446376204440720569753E+01,
		     er7 =-0.4957589496572501915214079952E+00,
		     er8 = 0.1664377182454986536961530415E+01,
		     er9 =-0.3503288487499736816886487290E+00,
		     er10= 0.3341791187130174790297318841E+00,
		     er11= 0.8192320648511571246570742613E-01,
		     er12=-0.2235530786388629525884427845E-01;
	double err=0.0,err2=0.0,sqr,sk,deno;

	// Calculate the contribution to the error from each variable
	for(int i=0;i<n;i++) {
		sk=1.0/(atoli+rtoli*max_d(fabs(w[i]),fabs(k5[i])));
		sqr=k4[i]-bhh1*k1[i]-bhh2*k9[i]-bhh3*k3[i];
		sqr*=sk;err2+=sqr*sqr;
		sqr=er1*k1[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+er10*k10[i]+er11*k2[i]+er12*k3[i];
		sqr*=sk;err+=sqr*sqr;
	}

	// Assemble the combination of the third-order and fifth-order
	// estimators
	deno=err+0.01*err2;
	return err*sqrt(1.0/(deno<=0.0?n:deno*n));
}

/** Detects whether the ODE system has become stiff.
 * \return True if the system has become stiff, false otherwise. */
bool dop853::detect_stiffness() {
	double sqr,stnum=0.0,stden=0.0;

	// Calculate the measure of stiffness, using previously computed RK
	// steps
	for(int i=0;i<n;i++) {
		sqr=k4[i]-k3[i];stnum+=sqr*sqr;
		sqr=k5[i]-ww1[i];stden+=sqr*sqr;
	}

	// If the stiffness criterion is met, and many recent steps have been
	// stiff, then bail out
	if(stden>0.0&&h*h*stnum>37.21*stden) {
		nonsti=0;
		iasti++;
		if(iasti==15) return true;
	}

	// If the system is not stiff after several tests, then reset the
	// stiffness counter
	if(++nonsti==6) iasti=0;
	return false;
}

/** Assuming that a regular Dormand--Prince step has been computed, this
 * function will carry out the necessary calculations of the seventh-order
 * dense output interpolation formulae. */
void dop853::dense_output() {
	int i;
	double ydiff,bspl;
	const double	c14 = 0.1E+00,
			c15 = 0.2E+00,
			c16 = 0.777777777777777777777777777778E+00,

			a141 =  5.61675022830479523392909219681E-2,
			a147 =  2.53500210216624811088794765333E-1,
			a148 = -2.46239037470802489917441475441E-1,
			a149 = -1.24191423263816360469010140626E-1,
			a1410 =  1.5329179827876569731206322685E-1,
			a1411 =  8.20105229563468988491666602057E-3,
			a1412 =  7.56789766054569976138603589584E-3,
			a1413 = -8.298E-3,

			a151 =  3.18346481635021405060768473261E-2,
			a156 =  2.83009096723667755288322961402E-2,
			a157 =  5.35419883074385676223797384372E-2,
			a158 = -5.49237485713909884646569340306E-2,
			a1511 = -1.08347328697249322858509316994E-4,
			a1512 =  3.82571090835658412954920192323E-4,
			a1513 = -3.40465008687404560802977114492E-4,
			a1514 =  1.41312443674632500278074618366E-1,

			a161 = -4.28896301583791923408573538692E-1,
			a166 = -4.69762141536116384314449447206E0,
			a167 =  7.68342119606259904184240953878E0,
			a168 =  4.06898981839711007970213554331E0,
			a169 =  3.56727187455281109270669543021E-1,
			a1613 = -1.39902416515901462129418009734E-3,
			a1614 =  2.9475147891527723389556272149E0,
			a1615 = -9.15095847217987001081870187138E0,

			d41 =-0.84289382761090128651353491142E+01,
			d46 = 0.56671495351937776962531783590E+00,
			d47 =-0.30689499459498916912797304727E+01,
			d48 = 0.23846676565120698287728149680E+01,
			d49 = 0.21170345824450282767155149946E+01,
			d410=-0.87139158377797299206789907490E+00,
			d411= 0.22404374302607882758541771650E+01,
			d412= 0.63157877876946881815570249290E+00,
			d413=-0.88990336451333310820698117400E-01,
			d414= 0.18148505520854727256656404962E+02,
			d415=-0.91946323924783554000451984436E+01,
			d416=-0.44360363875948939664310572000E+01,

			d51 = 0.10427508642579134603413151009E+02,
			d56 = 0.24228349177525818288430175319E+03,
			d57 = 0.16520045171727028198505394887E+03,
			d58 =-0.37454675472269020279518312152E+03,
			d59 =-0.22113666853125306036270938578E+02,
			d510= 0.77334326684722638389603898808E+01,
			d511=-0.30674084731089398182061213626E+02,
			d512=-0.93321305264302278729567221706E+01,
			d513= 0.15697238121770843886131091075E+02,
			d514=-0.31139403219565177677282850411E+02,
			d515=-0.93529243588444783865713862664E+01,
			d516= 0.35816841486394083752465898540E+02,

			d61= 0.19985053242002433820987653617E+02,
			d66=-0.38703730874935176555105901742E+03,
			d67=-0.18917813819516756882830838328E+03,
			d68= 0.52780815920542364900561016686E+03,
			d69=-0.11573902539959630126141871134E+02,
			d610= 0.68812326946963000169666922661E+01,
			d611=-0.10006050966910838403183860980E+01,
			d612= 0.77771377980534432092869265740E+00,
			d613=-0.27782057523535084065932004339E+01,
			d614=-0.60196695231264120758267380846E+02,
			d615= 0.84320405506677161018159903784E+02,
			d616= 0.11992291136182789328035130030E+02,

			d71 =-0.25693933462703749003312586129E+02,
			d76 =-0.15418974869023643374053993627E+03,
			d77 =-0.23152937917604549567536039109E+03,
			d78 = 0.35763911791061412378285349910E+03,
			d79 = 0.93405324183624310003907691704E+02,
			d710=-0.37458323136451633156875139351E+02,
			d711= 0.10409964950896230045147246184E+03,
			d712= 0.29840293426660503123344363579E+02,
			d713=-0.43533456590011143754432175058E+02,
			d714= 0.96324553959188282948394950600E+02,
			d715=-0.39177261675615439165231486172E+02,
			d716=-0.14972683625798562581422125276E+03;

	// Calculate the contributions to the dense output corrections using
	// the previously computed RK steps
	for(i=0;i<n;i++) {
		rc1[i]=w[i];
		ydiff=k5[i]-w[i];
		rc2[i]=ydiff;
		bspl=h*k1[i]-ydiff;
		rc3[i]=bspl;
		rc4[i]=ydiff-h*k4[i]-bspl;
		rc5[i]=d41*k1[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+d49*k9[i]+d410*k10[i]+d411*k2[i]+d412*k3[i];
		rc6[i]=d51*k1[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+d59*k9[i]+d510*k10[i]+d511*k2[i]+d512*k3[i];
		rc7[i]=d61*k1[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+d69*k9[i]+d610*k10[i]+d611*k2[i]+d612*k3[i];
		rc8[i]=d71*k1[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+d79*k9[i]+d710*k10[i]+d711*k2[i]+d712*k3[i];
	}

	// Carry out the next three function evaluations (steps 14 to 16)
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a141*k1[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]
			       +a1410*k10[i]+a1411*k2[i]+a1412*k3[i]+a1413*k4[i]);
	ff(t+c14*h,ww1,k10);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a151*k1[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]
			       +a1511*k2[i]+a1512*k3[i]+a1513*k4[i]+a1514*k10[i]);
	ff(t+c15*h,ww1,k2);
	for(i=0;i<n;i++) ww1[i]=w[i]+h*(a161*k1[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]
			       +a169*k9[i]+a1613*k4[i]+a1614*k10[i]+a1615*k2[i]);
	ff(t+c16*h,ww1,k3);
	nff+=3;

	// Use the newly computed steps to complete the calculation of the
	// dense output corrections
	for(i=0;i<n;i++) {
		rc5[i]=h*(rc5[i]+d413*k4[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
		rc6[i]=h*(rc6[i]+d513*k4[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
		rc7[i]=h*(rc7[i]+d613*k4[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
		rc8[i]=h*(rc8[i]+d713*k4[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
	}
}

/** Calculates a snapshot of the variables using the previously computed
 * seventh-order accurate dense output correction.
 * \param[in] ws an array of size n in which to store the variables.
 * \param[in] ti the time at which to compute the variables. */
void dop853::dense(double *ws,double ti) {
	double s=(ti-told)/h,s1=1.0-s;
	for(int i=0;i<n;i++) ws[i]=rc1[i]+s*(rc2[i]+s1*(rc3[i]+s*(rc4[i]
				  +s1*(rc5[i]+s*(rc6[i]+s1*(rc7[i]+s*rc8[i]))))));
}

/** Estimates the initial step size. It assumes that the k1 array has been set
 * to dw/dt.
 * \param[in] hmax the maximum step size allowable.
 * \param[in] posneg the direction of integration, set to +1 for forward
 *		     integration and -1 for backward integration.
 * \return The estimated step size. */
double dop853::hinit(double hmax,double posneg) {
	double dnf=0.,dny=0.,der2=0.,der12,h,h1,sk,sqr;
	int i;

	// Compute preliminary step size estimate
	for(i=0;i<n;i++) {
		sk=atoli+rtoli*fabs(w[i]);
		sqr=k1[i]/sk;
		dnf+=sqr*sqr;
		sqr=w[i]/sk;
		dny+=sqr*sqr;
	}
	h=min_d(dnf<=1E-10||dny<=1E-10?1.0E-6:sqrt(dny/dnf)*0.01,hmax)*posneg;

	// Perform an explicit Euler step
	for(i=0;i<n;i++) ww1[i]=w[i]+h*k1[i];
	ff(t+h,ww1,k2);nff++;

	// Estimate the second derivative of the solution
	for(i=0;i<n;i++) {
		sqr=(k2[i]-k1[i])/(atoli+rtoli*fabs(w[i]));
		der2+=sqr*sqr;
	}
	der2=sqrt(der2)/h;

	// The step size is computed such that h**8*max_d(norm(f0),norm(der2))=0.01
	der12=max_d(fabs(der2),sqrt(dnf));
	h1=der12<=1.0E-15?max_d(1.0E-6,fabs(h)*1.0E-3)
			 :pow(0.01/der12,0.125);
	return min_d(100.0*fabs(h),min_d(h1,hmax))*posneg;
}
