#include <iostream>

#include "margin-distrib.hpp"

using namespace MarginDistrib;
using std::cout;
using std::cerr;
using std::endl;


// ---------------------------------------------------------------------------
// class UnivarMargin

void UnivarMargin::inv_cdf(VectorI const & ranks, VectorD & cdfs)
{
	DimT N = ranks.size();
	if (N == 0)
		throw std::logic_error("called inv_cdf with zero-sized input vector!");
	if (cdfs.size() != N) {
		if (cdfs.size() > 0)
			cerr << "vector cdfs in inv_cdf has a wrong size -> resizing" << endl;
		cdfs.resize(N);
	}
	for (DimT i = 0; i < N; ++i) {
		cdfs(i) = inv_cdf(ranks(i), N);
	}
	if (fixEV) {
		fix_mean_std(cdfs, mean, (fixSD ? stDev : -1.0));
	}
}


// ---------------------------------------------------------------------------
// class MarginNormal

MarginNormal::MarginNormal(double const mu, double const sigma,
                           SamplePP const postP,
                           bool const condEVInv)
: UnivarMargin(postP, mu, sigma),
  invCdfF(mu, sigma), invCdf01F(0.0, 1.0), useCondEV(condEVInv),
  evMult(1.0 / sqrt(2 * 3.1415926535898) * sigma)
{}


double MarginNormal::inv_cdf(DimT const r, DimT const N) const
{
	if (useCondEV) {
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

double MarginSample::inv_cdf(double const p) const
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
