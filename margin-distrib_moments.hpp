/** \file
	definitions for using the HKW moment-matching heuristic.
**/

#ifndef MARGIN_DISTRIB_MOMENTS_HPP
#define MARGIN_DISTRIB_MOMENTS_HPP

#include "margins.hpp"

namespace MarginDistrib {

/// margin specified by four moments
/**
	\warning needs the moment-matching library by Michal Kaut
**/
class MarginMoments : public UnivarMargin {
private:
	VectorD moments;     ///< the first four moments
	int formOfMoms;      ///< format of the moments (as in the HKW code)
	VectorD sortedVals;  ///< the generated values, sorted

	/// generate and sort the scenario values
	void gen_scen(DimT const nSc);

	/// inverse CDF using ranks and number of scenarios
	/**
		Here we call the moment-matching method (since we do not know the
		number of scenarios on construction)

		\param[in] r  rank: 0,...,N-1
		\param[in] N  number of samples/scenarios
	**/
	boost::optional<double> inv_cdf_r(DimT const r, DimT const N) const override;

public:
	//MarginMoments(VectorD const & tgMoms, int const FoM = 0);
	MarginMoments(VectorD const & tgMoms, DimT const nSc, int const FoM = 0);

	/// constructor with parameters in a string stream, plus a number of scens.
	MarginMoments(std::istream & paramStr, DimT const nSc);
};

} // namespace


#endif // header guard
