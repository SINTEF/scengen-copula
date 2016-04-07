#include "margin-distrib.hpp"
#include "margin-distrib_moments.hpp"

#include <iostream>
#include <fstream>

using namespace MarginDistrib;
using std::cout;
using std::cerr;
using std::endl;
using std::string;


// ---------------------------------------------------------------------------
// generic routines

/*
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
*/

// ---------------------------------------------------------------------------
// class UnivarMargin

// creates the map - new derived classes have to be added here!
void UnivarMargin::init_name_map() {
#ifdef HAS_HKW
	add_to_map<MarginMoments>("moments");
	add_to_map<MarginMoments>("m");
#endif // HAS_HKW
#ifdef HAS_QUANTLIB
	add_to_map<MarginNormal>("normal");
	add_to_map<MarginNormal>("n");
#endif // HAS_QUANTLIB
	add_to_map<MarginTriang>("triang");
	add_to_map<MarginSample>("sample");
	add_to_map<MarginSample>("data");
	add_to_map<MarginStudent>("student");
	add_to_map<MarginLognormal>("lognormal");
	add_to_map<MarginExp>("exponential");
	add_to_map<MarginExp>("exp");
	add_to_map<MarginBeta>("beta");
	add_to_map<MarginPoisson>("poisson");
	add_to_map<MarginUniformC>("uniform");
}

// template for adding new entries to the-name-map
template <class T>
void UnivarMargin::add_to_map(string const & idStr)
{
	// throw an error if the element with the same id-string already exists
	if (nameMap.count(idStr) > 0) {
		throw std::logic_error("Tried to re-insert '" + idStr
		                       + "' into the name map.");
	}
	// convert to lower case
	string idStrLC = idStr;
	for (auto & ltr: idStrLC)
		ltr = std::tolower(ltr);
	// add to the map (using the new C++11 method)
	nameMap.emplace(idStrLC,
	                [](std::istream & paramStr, DimT const nSc) {
	                	return static_cast<UnivarMargin*>(new T (paramStr, nSc));
	                });
}

// creates a new object based on its name
/// \todo checks if shapeName is in the map
UnivarMargin * UnivarMargin::make_new(string const & margName,
                                      std::istream & paramStr,
                                      DimT const nSc)
{
	// create the map if it does not exist
	// this way, we do not need to create the map manually
	if (nameMap.size() == 0)
		init_name_map();
	// convert to lower case
	string margNameLC = margName;
	for (auto & ltr: margNameLC)
		ltr = std::tolower(ltr);
	if (nameMap.count(margNameLC) == 0) {
		string errMsg = "Unknown margin-distribution type '" + margName + "'. "
		                + "Supported distribution ids are:";
		for (auto nm : nameMap)
			errMsg += " " + nm.first;
		throw std::runtime_error (errMsg + ".");
	}
	return UnivarMargin::nameMap.at(margNameLC)(paramStr, nSc);
}

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

// initialize the map - required!
UnivMargNameMapT UnivarMargin::nameMap;

#ifdef HAS_QUANTLIB
// ---------------------------------------------------------------------------
// class MarginNormal

MarginNormal::MarginNormal(double const mu, double const sigma,
                           SamplePP const postP, bool const condEVInv)
: UnivarMargin(postP, mu, sigma),
  invCdfF(mu, sigma), invCdf01F(0.0, 1.0), useCondEV(condEVInv),
  evMult(sigma / sqrt(2 * 3.1415926535898))
{}

MarginNormal::MarginNormal(std::istream & paramStr, DimT const nSc,
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
#endif // HAS_QUANTLIB


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

// constructor with parameters in a string stream
/// \todo Take advantage of nSc, or remove it from the parameter list
MarginSample::MarginSample(std::istream & paramStr, DimT const nSc,
                           SamplePP const postP)
: UnivarMargin(postP)
{
	std::string inFName;
	DimT inFCol;

	paramStr >> inFName;
	if (paramStr.fail()) {
		throw std::runtime_error
			("error while reading parameters for a sample margin");
	}
	std::ifstream inFStr(inFName.c_str());
	if (!inFStr) {
		throw std::ios_base::failure("Could not open input file `" + inFName + "'!");
	}
	paramStr >> inFCol;
	if (paramStr.fail()) {
		// no second parameter -> data in a vector format
		inFStr >> sortedS; // assumes dimension given in the file
	} else {
		// second parameter is the column index -> read as a matrix
		MatrixD tmpMat;
		inFStr >> tmpMat; // assumes dimensions given in the file
		sortedS = ublas::column(tmpMat, inFCol-1); // assume counting from 1
	}
	if (!inFStr.good()) {
		throw std::ios_base::failure("Error while reading file `" + inFName + "'!");
	}
	inFStr.close();
	nPts = sortedS.size();

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

MarginTriang::MarginTriang(std::istream & paramStr, DimT const nSc, bool const useExtremes)
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

MarginExp::MarginExp(std::istream & paramStr, DimT const nSc)
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

MarginBeta::MarginBeta(std::istream & paramStr, DimT const nSc)
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

MarginLognormal::MarginLognormal(std::istream & paramStr, DimT const nSc)
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

MarginStudent::MarginStudent(std::istream & paramStr, DimT const nSc)
: UnivarMargin(), p2Dist()
{
	paramStr >> dof;
	p2Dist.reset(new boost::math::students_t_distribution<>(dof));
	if (paramStr.good()) {
		paramStr >> mean >> stDev;
		if (! paramStr.fail()) {
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

MarginPoisson::MarginPoisson(std::istream & paramStr, DimT const nSc)
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
MarginUniformC::MarginUniformC(std::istream & paramStr, DimT const nSc)
: UnivarMargin()
{
	paramStr >> a >> b;
	supLen = b - a;
}

boost::optional<double> MarginUniformC::inv_cdf(double const p) const
{
	return a + p * supLen;
}
