#ifndef MARGIN_DISTRIB_HPP
#define MARGIN_DISTRIB_HPP

#include <ql/math/distributions/normaldistribution.hpp> // for normal distrib.

#include "common.hpp"

namespace MarginDistrib {

// ----------------------------------------------------------------------------
class UnivarMargin {
protected:
	bool fixEV;   ///< should we fix the exp. val. (when used on vector)?
	bool fixSD;   ///< should we fix the std. dev. (when used on vector)?
	double mean;  ///< mean
	double stDev; ///< standard deviation

public:
	/// postprocess applied to inv_cdf of a whole margin
	enum SamplePP { PP_None, PP_fixMean, PP_fixBoth };

	UnivarMargin(SamplePP const postP = PP_None,
	             double const EV = 0.0, double const SD = -1.0)
	: fixEV(postP != PP_None), fixSD(postP == PP_fixBoth), mean(EV), stDev(SD)
	{}
	virtual ~UnivarMargin() {}

	virtual double inv_cdf(double const p) const = 0;

	virtual double inv_cdf(DimT const r, DimT const N) const {
		return inv_cdf((static_cast<double>(r) + 0.5) / N);
	}

	virtual void inv_cdf(VectorI const & ranks, VectorD & cdfs);

	typedef boost::shared_ptr<UnivarMargin> Ptr;
};


// ----------------------------------------------------------------------------
class MarginNormal : public UnivarMargin {
private:
	/// object that computes the inverse - mean and std. dev. given at constr.
	QuantLib::InverseCumulativeNormal invCdfF;

	/// \name inverse using conditional means
	/**
		If we have a rank \c r from \c 0,...,\c N-1, it represents interval
		\c (r/N, r+1/N). The default approach is to return inverse-cdf of the
		middle of the interval, i.e. \c (r+0.5)/N. \n
		One possible alternative is to compute the expected value of the variable,
		conditional on the fact that its cdf is in the interval, i.e. that
		the variable is in \f$ \bigl( F^{-1}(r/N), F^{-1}((r+1)/N) \bigr) \f$.
		This is a bit more complicated, but guarantees that the whole sample
		will have the correct expected value. \n
		In this sense, the default method corresponds to the conditional \em
		median of the same interval.
	**/
	///@{
		QuantLib::InverseCumulativeNormal invCdf01F; ///< N(0,1) inverse - needed
		bool useCondEV;      ///< use conditional means?
		double const evMult; ///< constant used in the calculation
	///@}

public:
	MarginNormal(double const mu, double const sigma,
	             SamplePP const postP = PP_None,
	             bool const condEVInv = false);
	~MarginNormal() {}

	double inv_cdf(double const p) const { return invCdfF(p); }

	double inv_cdf(DimT const r, DimT const N) const;
};


// ----------------------------------------------------------------------------
class MarginSample : public UnivarMargin {
private:
	VectorD sortedS; ///< sorted sample
	DimT nPts;       ///< number of sample points

public:
	MarginSample(VectorD const & sample, SamplePP const postP = PP_None);

	double inv_cdf(double const p) const;
};


// ----------------------------------------------------------------------------
/// triangular distribution
class MarginTriang : public UnivarMargin {
private:
	double min;
	double max;
	double mode;
	bool useMinMax; ///< should scenarios always include min and max?

	double modeCdf; ///< used for inverse cdf
	double valMin;  ///< used for inverse cdf
	double valMax;  ///< used for inverse cdf

public:
	MarginTriang(double const minV, double const maxV, double const modeV,
	             bool const useExtremes = false);

	double inv_cdf(double const p) const;

	double inv_cdf(DimT const r, DimT const N) const;
};


} // namespace MarginDistrib

#endif
