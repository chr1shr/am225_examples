#ifndef TGMGPP_PREDICT_HH
#define TGMGPP_PREDICT_HH

#include "tgmg_config.hh"

struct tgmg_predict {
	const int mult;
	const int decay;
	const int max_thresh;
	const int extra_iters;
	/* A current estimate of the number of V-cycles required for
	 * convergence, multiplied by the mgs_mult constant. */
	int lim;
	/** A running counter of the number of multigrid solves that have been
	 * performed, used for printing diagnostic information about the
	 * avergage number of V-cycles. */
	int solves;
	/** A running counter of the number of multigrid V-cycles that have
	 * been performed. */
	int vcount;
	tgmg_predict() : mult(tgmg_predict_mult), decay(tgmg_predict_decay),
		max_thresh(mult*tgmg_predict_max_iters), extra_iters(tgmg_predict_extra_iters),
		lim(mult*tgmg_predict_init), solves(0), vcount(0) {}
	inline void add_iters(int n) {
		solves++;
		vcount+=n;
	}
	inline float avg_iters() {
		float ans=static_cast<float>(vcount)/static_cast<float>(solves);
		solves=vcount=0;
		return ans;
	}
};

#endif
