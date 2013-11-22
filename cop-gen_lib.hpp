#ifndef COPGEN_LIB_HPP
#define COPGEN_LIB_HPP

/// \file cop-gen_lib.hpp
/**
	This file provides the interface for the cop-gen library
**/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ScenGenCop {

//namespace ublas = boost::numeric::ublas; // shortcut name

int gen_scen_normal(boost::numeric::ublas::vector<double> const & tgMean,
                    boost::numeric::ublas::vector<double> const & tgStD,
                    boost::numeric::ublas::matrix<double> const & tgCorr,
                    unsigned const nSc,
                    boost::numeric::ublas::matrix<double> & scens);

} // namespace ScenGenCop

#endif // header guard
