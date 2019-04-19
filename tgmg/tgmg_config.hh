#ifndef TGMGPP_CONFIG_HH
#define TGMGPP_CONFIG_HH

/** The maximum number of grid levels to allocate storage space for. */
const int tgmg_max_levels=16;

const int jacobi_per_loop=100;
const int max_jacobi_loops=400;

const int gs_per_loop=20;
const int max_gs_loops=2000;

const int sor_per_loop=20;
const int max_sor_loops=1000;

const int v_cycle_per_loop=4;
const int max_v_cycle_loops=1000;

const int fmg_per_loop=2;
const int max_fmg_loops=1000;

/** The number of grid points at which to stop introducing coarser grids. */
const int tgmg_grid_min=10;

/** The default initial number of V-cycles to use when solving the multigrid problem. */
const int tgmg_predict_init=8;

/** A default multiplier to apply to the internal counter for tracking the
 * number of V-cycles required. After each solve, the counter is decreased by
 * one (plus an additional multiplicative decay factor). The smaller the
 * multiplier, the more aggressive the solver is at trying to reduce the
 * required number of V-cycles. */
const int tgmg_predict_mult=16;

/** A default decay factor to try and lower the number of multigrid cycles
 * performed at each time. */
const int tgmg_predict_decay=32;

/** The default maximum number of V-cycles to perform before bailing out with an error.
 */
const int tgmg_predict_max_iters=2048;

/** The default number of V-cycles to perform after reaching the tolerance, to
 * achieve further accuracy. */
const int tgmg_predict_extra_iters=1;

#ifdef _OPENMP
#include "omp.h"
const double tgmg_time_group=omp_get_wtick()*1000;
const double tgmg_time_max=omp_get_wtick()*3200;
const double tgmg_f_bail_factor=3;
#endif

const double tgmg_retune_interval=3600.;
const double tgmg_default_resource_factor=0.1;

/** A status value to return if a fatal error is encountered. */
#define TGMGPP_ERROR 1

#endif
