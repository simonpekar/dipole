#include <TMath.h>

/* constants */
const Double_t zeta_3 = 1.2020569031595942854;

/* general parameters */
const Double_t delta = 0.001; /* UV cutoff (between 10^-6 and 10^-2) */
const Double_t y_max = 3.; /* maximum rapidity in evolution (typ 3.) */
const UInt_t N = 10000; /* number of independant events (typ 10000) */

/* IR cutoff parameters */
const Double_t R = 2.; /* IR cutoff (typ 2.), 0. for no cutoff */
const UInt_t IR_type = 1; /* type of IR cutoff */

/* fluctuations calculation parameters */
const Bool_t compute_fluctuations = true; /* compute fluctutations */
const Double_t r_s = exp(-5); /* minimum size to consider (~ delta) */

/* common ancestors parameters */
const Bool_t compute_ancestors = false; /* compute common ancestors ? */
const UInt_t k_leaves = 1; /* compute the common ancestor of k random leaves (typ 2-10) */
const UInt_t m_factor = 6; /* consider only events with a multiplicity greater than mfactor * mean multiplicity (typ 2) */
