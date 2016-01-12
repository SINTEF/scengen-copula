/// declaration of 2D copulas that depend on QuantLib library
/**
	This makes it easier to create versions of the code that does neither
	include or depend on the QuantLib library

	For the moment, this file includes normal and student copulas
**/

#ifndef COP_2D_INFO_QUANTLIB_HPP
#define COP_2D_INFO_QUANTLIB_HPP

#include "cop2Dinfo.hpp"

// QuantLib libraries used for the normal copula
#include "external/QuantLib/math/distributions/bivariatenormaldistribution.hpp"
#include "external/QuantLib/math/distributions/bivariatestudenttdistribution.hpp"
#include "external/QuantLib/math/distributions/studenttdistribution.hpp"

//#include <boost/optional.hpp>
//#include <iostream>
//#include <map>

namespace Copula2D {

/// Normal copula
/**
	\note I have not found any library providing the copula cdf, so I use the
	      bivariate distribution and then inverse of the marginal cdfs.
	      For the former, the only implementation with suitable license is from
	      QuantLib; for the latter, we use QuantLib, but boost have the function
	      as well...
**/
class Cop2DNormal : public Cop2DComp {
private:
	double correl;  ///< correlation

	/// \name QuantLib objects that compute normal CDFs etc.
	///@{
		/// inverse cdf (by default standard normal distrib.)
		QuantLib::InverseCumulativeNormal N01InvCdf;
		/// bivariate cdf (standardized, needs correlation)
		/**
			\note Must use a pointer, since this needs correlation on construction
			      and not all class constructors have this as a parameter.
		**/
		std::unique_ptr<QuantLib::BivariateCumulativeNormalDistribution> p2N01Cdf2D;
	///@}

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DNormal(double const rho);

	Cop2DNormal(std::istream & paramStr);
};


/// Bivariate Student-t copula
/**
	I have not found any library (with suitable license) providing bivariate
	student copula, nor bivariate student distribution. As a result, I had
	implemented the distribution myself and then compute the copula cdf()
	using the univariate inverse CDF from QuantLib.
	(There is also one in boost, but we are already using QuantLib for other
	things...)
**/
class Cop2DStudent : public Cop2DComp {
private:
	double correl;  ///< correlation
	unsigned dof;   ///< degree of freedom

	/// \name QuantLib objects that compute CDFs etc.
	/**
		\note Must use pointers, since the object needs dof on construction
		      and not all class constructors have this as a parameter.
	**/
	///@{
		///class providing bivariate student cdf
		std::unique_ptr<QuantLib::BivariateCumulativeStudentDistribution> p2tCdf2D;

		/// QuantLib object for the inverse cdf
		std::unique_ptr<QuantLib::InverseCumulativeStudent> p2tInvCdf;
	///@}

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DStudent(unsigned degF, double const rho);

	Cop2DStudent(std::istream & paramStr);
};

} // namespace Copula2D

#endif // header guard
