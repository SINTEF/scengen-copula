#ifndef MARGIN_DISTRIB_HPP
#define MARGIN_DISTRIB_HPP

#include "common.hpp"

#include <ql/math/distributions/normaldistribution.hpp> // for normal distrib.

#include <map>
#include <boost/optional.hpp>


/// classes and methods related to the marginal distributions
namespace MarginDistrib {

/// \name objects for the margin-type name map, used in the main code
///@{
	/// enum for the known multivariate margin specification types
	enum class MargDistribID {moments, normal, sample, triang, triangX, unknown};

	/// type for the copula map
	typedef std::map<std::string, MargDistribID> DistribNameMapT;

	/// this fills the copula map with the known copula types
	void make_distrib_name_map(DistribNameMapT & dMap);
///@}


// ----------------------------------------------------------------------------
/// base class for all the marginal-distribution classes (pure virtual)
class UnivarMargin {
protected:
	bool fixEV;   ///< should we fix the exp. val. (when used on vector)?
	bool fixSD;   ///< should we fix the std. dev. (when used on vector)?
	double mean;  ///< mean
	double stDev; ///< standard deviation

	/// inverse CDF using ranks and number of scenarios
	/**
		This enables the margin-distribution classes to implement inverse CDF
		using ranks, should this be easier/faster than using values from (0,1).
		By default, this returns an empty value.

		\param[in] r  rank: 0,...,N-1
		\param[in] N  number of samples/scenarios
	**/
	virtual boost::optional<double> inv_cdf_r(DimT const r, DimT const N) const {
		return boost::optional<double>();
	}

	/// inverse CDF with single argument (the percentile)
	/**
		This should be overwritten by classes where this version is preferable.
		By default, this returns an empty value.

		\param[in] p the percentile
	**/
	virtual boost::optional<double> inv_cdf(double const p) const {
		return boost::optional<double>();
	}

public:
	/// postprocess applied to inv_cdf of a whole margin
	enum SamplePP { PP_None, PP_fixMean, PP_fixBoth };

	UnivarMargin(SamplePP const postP = PP_None,
	             double const EV = 0.0, double const SD = -1.0)
	: fixEV(postP != PP_None), fixSD(postP == PP_fixBoth), mean(EV), stDev(SD)
	{}
	virtual ~UnivarMargin() {}

	/// inverse CDF
	/**
		By default, this is a wrapper around the protected methods.
		It uses the rank-based version if it is available,
		then the percentile-based one. (It will throw if none is available.)

		\param[in] r  rank: 0,...,N-1
		\param[in] N  number of samples/scenarios
	**/
	virtual double inv_cdf(DimT const r, DimT const N) const;

	/// inverse CDF, for all scenarios at once
	/**
		The default implementation simply calls the univariate version
		element-by-element.

		\param[in]   ranks  input values (percentiles)
		\param[out] values  output values (actual values from the distribution)
	**/
	virtual void inv_cdf(VectorI const & ranks, VectorD & values);

	typedef boost::shared_ptr<UnivarMargin> Ptr;
};


// ----------------------------------------------------------------------------
/// margin with normal distribution
/**
	uses objects from the open-source QuantLib library
**/
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

	boost::optional<double> inv_cdf(double const p) const override
	{ return invCdfF(p); }

public:
	MarginNormal(double const mu, double const sigma,
	             SamplePP const postP = PP_None,
	             bool const condEVInv = false);
	~MarginNormal() {}

	double inv_cdf(DimT const r, DimT const N) const override;
};


// ----------------------------------------------------------------------------
/// margin specified by a sample, i.e. a list of equiprobable values
/**
	The inverse CDF uses a linear interpolation of the inverse empirical CDF
**/
class MarginSample : public UnivarMargin {
private:
	VectorD sortedS; ///< sorted sample
	DimT nPts;       ///< number of sample points

	/// inverse CDF using ranks and number of scenarios
	/**
		This is used in the case the number of scenarios is the same as the
		number of samples points \a nPts.

		\param[in] r  rank: 0,...,N-1
		\param[in] N  number of samples/scenarios
	**/
	boost::optional<double> inv_cdf_r(DimT const r, DimT const N) const override
	{ return (N == nPts ? sortedS[r] : boost::optional<double>()); }

	boost::optional<double> inv_cdf(double const p) const override;

public:
	MarginSample(VectorD const & sample, SamplePP const postP = PP_None);
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

	void guts_of_constructor();

	boost::optional<double> inv_cdf(double const p) const override;

public:
	MarginTriang(double const minV, double const maxV, double const modeV,
	             bool const useExtremes = false);

	MarginTriang(std::stringstream & paramStr, bool const useExtremes = false);

	double inv_cdf(DimT const r, DimT const N) const override;
};


} // namespace MarginDistrib

#endif
