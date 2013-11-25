#ifndef COPGEN_LIB_HPP
#define COPGEN_LIB_HPP

/// \file cop-gen_lib.hpp
/**
	This file provides the interface for the cop-gen library
**/

#include "dll_export_def.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace ScenGenCop {

#ifdef NDEBUG
	unsigned defVerb = 3; // default verbosity of release code [0..10]
#else
	unsigned defVerb = 7; // default verbosity of debug code [0..10]
#endif


/// generate scenarios from normal distribution
/**
	\param[in] tgMean  vector of target means
	\param[in]  tgStD  vector of target standard deviations
	\param[in] tgCorr  matrix of target correlations
	\param[in]    nSc  number of scenarios to generate (equiprobable)
	\param[out] scens  matrix where we put the results
	\param[in] varInC  if true, output is [nSc * nVar], otherwise [nVar * nSc]
	\param[in]   verb  verbosity level (output); must be from {0,..,10}

	\return  the distance from the copula generator
**/
DLL_PUBLIC
double gen_scen_normal(boost::numeric::ublas::vector<double> const & tgMean,
                       boost::numeric::ublas::vector<double> const & tgStD,
                       boost::numeric::ublas::matrix<double> const & tgCorr,
                       unsigned const nSc,
                       boost::numeric::ublas::matrix<double> & scens,
                       bool const varInC = true,
                       unsigned const verb = defVerb);

} // namespace ScenGenCop

#endif // header guard
