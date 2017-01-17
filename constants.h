/* constants */
const Double_t zeta_3 = 1.2020569031595942854;

/* general parameters */
const Double_t delta = 0.001; /* UV cutoff (between 10^-6 and 10^-2) */
const Double_t y_max = 3.; /* maximum rapidity in evolution (typ 3.) */
const UInt_t N = 1000; /* number of independant events (typ 10000) */
const Double_t precision = 1./sqrt(N); /* for integration */

/* IR cutoff parameters */
const Bool_t IR = false; /* activate IR cutoff */
const Double_t R = 2.; /* IR cutoff for the Gaussian form (typ 2.) */

/* fluctuations calculation parameters */
const Bool_t compute_fluctuations = true; /* compute fluctutations */
const Double_t r_s = 0.005; /* minimum size to consider (~ delta) */

/* common ancestors parameters */
const Bool_t compute_ancestors = true; /* compute common ancestors ? */
const UInt_t k_leaves = 2; /* compute the common ancestor of k random leaves (typ 2-10) */
const UInt_t m_factor = 5; /* consider only events with a multiplicity greater than mfactor * mean multiplicity (typ 2) */
