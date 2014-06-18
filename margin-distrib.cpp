#include <iostream>

#include "margin-distrib.hpp"

using namespace MarginDistrib;
using std::cout;
using std::cerr;
using std::endl;


// ---------------------------------------------------------------------------
// generic routines

void MarginDistrib::make_distrib_name_map(DistribNameMapT & dMap) {
	// note: when listing the map, it goes from last to first
	dMap["m"] = MargDistribID::moments;
	dMap["moments"] = MargDistribID::moments;
	dMap["sample"] = MargDistribID::sample;
	dMap["n"] = MargDistribID::normal;
	dMap["normal"] = MargDistribID::normal;
	dMap["t"] = MargDistribID::triang;
	dMap["triang"] = MargDistribID::triang;
	dMap["tX"] = MargDistribID::triangX;
	dMap["triangX"] = MargDistribID::triangX;
	dMap["exponential"] = MargDistribID::exponential;
	dMap["exp"] = MargDistribID::exponential;
	dMap["beta"] = MargDistribID::beta;
	dMap["lognormal"] = MargDistribID::lognormal;
	dMap["poisson"] = MargDistribID::poisson;
	dMap["Poisson"] = MargDistribID::poisson;
	dMap["uniform"] = MargDistribID::uniform;
}

// ---------------------------------------------------------------------------
// class UnivarMargin

// inverse CDF
double UnivarMargin::inv_cdf(DimT const r, DimT const N) const
{
	boost::optional<double> x = inv_cdf_r(r, N);
	if (! x) {
		// inv_cdf_r() did returned empty value
		x = inv_cdf((static_cast<double>(r) + 0.5) / N);
		if (! x)
			throw std::logic_error("no inverse cdf defined for the margin!");
	}
	return x.get();
}

void UnivarMargin::inv_cdf(VectorI const & ranks, VectorD & values)
{
	DimT N = ranks.size();
	if (N == 0)
		throw std::logic_error("called inv_cdf with zero-sized input vector!");
	if (values.size() != N) {
		if (values.size() > 0)
			cerr << "vector values in inv_cdf has a wrong size -> resizing" << endl;
		values.resize(N);
	}
	for (DimT i = 0; i < N; ++i) {
		values(i) = inv_cdf(ranks(i), N);
	}
	if (fixEV) {
		fix_mean_std(values, mean, (fixSD ? stDev : -1.0));
	}
}


// ---------------------------------------------------------------------------
// class MarginNormal

MarginNormal::MarginNormal(double const mu, double const sigma,
                           SamplePP const postP, bool const condEVInv)
: UnivarMargin(postP, mu, sigma),
  invCdfF(mu, sigma), invCdf01F(0.0, 1.0), useCondEV(condEVInv),
  evMult(sigma / sqrt(2 * 3.1415926535898))
{}

MarginNormal::MarginNormal(std::stringstream & paramStr,
                           SamplePP const postP, bool const condEVInv)
: UnivarMargin(postP, 0, 1),
  invCdfF(0, 1), invCdf01F(0.0, 1.0), useCondEV(condEVInv),
  evMult(1.0 / sqrt(2 * 3.1415926535898))
{
	double mu, sigma;
	paramStr >> mu >> sigma;
	invCdfF = QuantLib::InverseCumulativeNormal(mu, sigma);
	evMult *= sigma;
}


double MarginNormal::inv_cdf(DimT const r, DimT const N) const
{
	if (useCondEV) {
		// For normal distribution, the conditional mean on interval [a,b] is
		// given by: -(f(b) - f(a)) / (F(b) - F(b)) = (f(a) - f(b)) / Pr([a,b]).
		// This is a special property of normal distribution, following from the
		// fact that integral of $x e^(-x^2/2)$ is $-e^(-x^2/2)$.
		double p1 = static_cast<double>(r) / N;
		double p2 = static_cast<double>(r + 1) / N;
		double F1 = (r == 0     ? 0.0 : exp(-0.5 * pow(invCdf01F(p1), 2)));
		double F2 = (r == N - 1 ? 0.0 : exp(-0.5 * pow(invCdf01F(p2), 2)));
		return N * evMult * (F1 - F2) + mean;
	} else {
		return invCdfF((static_cast<double>(r) + 0.5) / N);
	}
}


// ---------------------------------------------------------------------------
// class MarginSample

MarginSample::MarginSample(VectorD const & sample, SamplePP const postP)
: UnivarMargin(postP), sortedS(sample), nPts(sample.size())
{
	if (postP != PP_None) {
		mean = vec_mean(sortedS);
		if (postP == PP_fixBoth) {
			stDev = vec_std_dev(sortedS, mean); // using MLE formulas
		}
	}

	std::sort(sortedS.begin(), sortedS.end());
}

boost::optional<double> MarginSample::inv_cdf(double const p) const
{
	// the sample values are at points (i + 0.5)/nPts, for 0 <= i < nPts
	DimT i0;
	if (p < 0.5 / nPts) {
		i0 = 0;
	} else {
		i0 = floor(p * nPts - 0.5);
	}
	if (i0 == nPts - 1)
		i0--;
	double p0 = (i0 + 0.5) / nPts;
	double pos = (p - p0) * nPts; // 1/(p1 - p0) = 1/(1 / nPts) = nPts
	return pos * sortedS(i0) + (1 - pos) * sortedS(i0 + 1);
}


// ---------------------------------------------------------------------------
// class MarginTriang

void MarginTriang::guts_of_constructor()
{
	modeCdf = (mode - min) / (max - min);
	valMin = (max - min) * (mode - min);
	valMax = (max - min) * (max - mode);

	TRACE(TrDetail, "Created a new triangular margin with parameters ("
	                << min << ", " << max << ", " << mode << ")");
}

MarginTriang::MarginTriang(double const minV, double const maxV,
                           double const modeV, bool const useExtremes)
: UnivarMargin(), min(minV), max(maxV), mode(modeV), useMinMax(useExtremes)
{
	guts_of_constructor();
}

MarginTriang::MarginTriang(std::stringstream & paramStr, bool const useExtremes)
: UnivarMargin(), useMinMax(useExtremes)
{
	paramStr >> min >> max >> mode;
	guts_of_constructor();
}

boost::optional<double> MarginTriang::inv_cdf(double const p) const
{
	double x = (p < modeCdf ? min + sqrt(p * valMin)
	                        : max - sqrt((1-p) * valMax));
	return x;
}

double MarginTriang::inv_cdf(DimT const r, DimT const N) const
{
	double p = (useMinMax ? (static_cast<double>(r)) / (N - 1)
	                      : (static_cast<double>(r) + 0.5) / N);
	return inv_cdf(p).get();
}


// ---------------------------------------------------------------------------
// class MarginExp

MarginExp::MarginExp(double const rate)
: UnivarMargin(), lambda(rate)
{
	if (lambda <=0)
		throw std::range_error("parameter lambda is out of range");
}

MarginExp::MarginExp(std::stringstream & paramStr)
: UnivarMargin()
{
	paramStr >> lambda;
	if (lambda <=0)
		throw std::range_error("parameter lambda is out of range");
}

boost::optional<double> MarginExp::inv_cdf(double const p) const
{
	return -log(1 - p) / lambda;
}


// ---------------------------------------------------------------------------
// class MarginBeta

// constructor for the standard beta distribution
MarginBeta::MarginBeta(double const pAlpha, double const pBeta)
: UnivarMargin(), alpha(pAlpha), beta(pBeta), scaled(false),
  p2Dist(new boost::math::beta_distribution<>(pAlpha, pBeta))
{}

/// constructor for the standard beta distribution
MarginBeta::MarginBeta(double const pAlpha, double const pBeta,
                       double const min, double const max)
: UnivarMargin(), alpha(pAlpha), beta(pBeta), scaled(true),
  a(min), b(max), supLen(max - min),
  p2Dist(new boost::math::beta_distribution<>(pAlpha, pBeta))
{}

MarginBeta::MarginBeta(std::stringstream & paramStr)
: UnivarMargin(), p2Dist()
{
	paramStr >> alpha >> beta;
	p2Dist.reset(new boost::math::beta_distribution<>(alpha, beta));
	if (!paramStr.eof()) {
		scaled = true;
		paramStr >> a >> b;
		supLen = b - a;
	}
}

boost::optional<double> MarginBeta::inv_cdf(double const p) const
{
	double x = quantile(*p2Dist, p);
	if (scaled)
		x = a + supLen * x;
	return x;
}


// ---------------------------------------------------------------------------
// class MarginLognormal

// constructor for the standard Lognormal distribution
MarginLognormal::MarginLognormal(double const loc, double const scale)
: UnivarMargin(), mu(loc), sigma(scale),
  p2Dist(new boost::math::lognormal_distribution<>(loc, scale))
{}

MarginLognormal::MarginLognormal(std::stringstream & paramStr)
: UnivarMargin(), p2Dist()
{
	paramStr >> mu >> sigma;
	p2Dist.reset(new boost::math::lognormal_distribution<>(mu, sigma));
}

boost::optional<double> MarginLognormal::inv_cdf(double const p) const
{
	return quantile(*p2Dist, p);
}


// ---------------------------------------------------------------------------
// class MarginStudent

// constructor for the standard Lognormal distribution
MarginStudent::MarginStudent(unsigned const v, SamplePP const postP)
: UnivarMargin(postP),
  p2Dist(new boost::math::students_t_distribution<>(v))
{}

MarginStudent::MarginStudent(unsigned const v,
                             double const mu, double const sigma)
: UnivarMargin(PP_fixBoth, mu, sigma), scaled(true),
  p2Dist(new boost::math::students_t_distribution<>(v))
{}

MarginStudent::MarginStudent(std::stringstream & paramStr)
: UnivarMargin(), p2Dist()
{
	paramStr >> dof;
	p2Dist.reset(new boost::math::students_t_distribution<>(dof));
	if (!paramStr.eof()) {
		paramStr >> mean >> stDev;
		if (paramStr.good()) {
			scaled = true;
		}
	}
}

boost::optional<double> MarginStudent::inv_cdf(double const p) const
{
	return quantile(*p2Dist, p);
}


// ---------------------------------------------------------------------------
// class MarginPoisson

// constructor for the standard Poisson distribution
MarginPoisson::MarginPoisson(double const mean)
: UnivarMargin(), lambda(mean),
  p2Dist(new boostPoissonDistrib(mean))
{}

MarginPoisson::MarginPoisson(std::stringstream & paramStr)
: UnivarMargin(), p2Dist()
{
	paramStr >> lambda;
	p2Dist.reset(new boostPoissonDistrib(lambda));
}

boost::optional<double> MarginPoisson::inv_cdf(double const p) const
{
	return quantile(*p2Dist, p);
}


// ---------------------------------------------------------------------------
// class MarginUniformC

// constructor for the standard beta distribution
MarginUniformC::MarginUniformC(double const min, double const max)
: UnivarMargin(), a(min), b(max), supLen(max - min)
{}

// constructor with parameters in a string stream
MarginUniformC::MarginUniformC(std::stringstream & paramStr)
: UnivarMargin()
{
	paramStr >> a >> b;
	supLen = b - a;
}

boost::optional<double> MarginUniformC::inv_cdf(double const p) const
{
	return a + p * supLen;
}
