// These functions implement d0 normalization. The d0 for final TM-score
// output is implemented by parameter_set4final. For both RNA alignment
// and protein alignment, using d0 set by parameter_set4search yields
// slightly better results during initial alignment-superposition iteration.

#pragma once
#include "io.h"

void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0);

void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search);

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type);

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search);