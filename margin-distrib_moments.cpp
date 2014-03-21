#include "margin-distrib_moments.hpp"

/// namespace for the external code for moment-matching scenario generation
/**
	The code implements heuristic from paper
	'<em>A Heuristic for Moment-matching Scenario Generation</em>'
	by Kjetil Høyland, Michal Kaut and Stein W. Wallace, published in
	<em>Computational Optimization and Applications</em>, 24 (2-3), pp. 169–185,
	2003; <a href="http://dx.doi.org/doi:10.1023/A:1021853807313">
	doi:10.1023/A:1021853807313</a>. It is distributed under the
	<a href="http://www.eclipse.org/legal/epl-v10.html">
	Eclipse Public License</a>.
**/
namespace ExtScenGenHKW {
	extern "C" {
		#define NO_DLL_DEFS // not using DLLs -> ignore DLL export/import macros
		#include "external/scen-gen_HKW/HKW_sg.h"
	}
}

#include <fstream>
#include <ios>
#include <algorithm>

using namespace MarginDistrib;
using std::cout;
using std::endl;
using std::cerr;

// --------------------------------------------------------------------------
// class MarginMoments

MarginMoments::MarginMoments(VectorD const & tgMoms, DimT const nSc,
                             int const FoM)
: UnivarMargin(), moments(tgMoms), formOfMoms(FoM)
{
	gen_scen(nSc);
}


boost::optional<double> MarginMoments::inv_cdf_r(DimT const r, DimT const N) const
{
	if (sortedVals.size() != N) {
		if (sortedVals.size() == 0) {
			MSG (TrError,
				"scenario values must be generated before calling inv_cdf_r()!");
		} else {
			MSG (TrError,
				"generated scenario values have wrong size in inv_cdf_r()!");
		}
		exit(EXIT_FAILURE);
	}
	return sortedVals(r);
}


void MarginMoments::gen_scen(DimT const nSc)
{
	if (sortedVals.size() == nSc) {
		MSG (TrWarn, "re-generating scenarios - this should not happen");
	} else {
		if (sortedVals.size() > nSc)
			MSG (TrWarn, "re-generating scenarios with different size - strange!");
		sortedVals.resize(nSc, false);
	}

	//! scengen_HKW() uses either its own matrix classes, or C pointers
	//! -> have to create temp. objects with the right format
	DimT i;
	double * hkwTgMoms[4];   ///< target moments
	for (i = 0; i < 4; ++i) {
		hkwTgMoms[i] = &moments(i);
	}
	double * hkwTgCorrs[1];  ///< target correlations
	double dummyCorr = 1.0;
	hkwTgCorrs[0] = &dummyCorr;
	double * prob = nullptr;  ///< target probabilities (NULL for equiprob.)
	double * hkwOutSc[1];     ///< the output array
	hkwOutSc[0] = &sortedVals(0); // write directly to the vector - check!
	int hkwOutLvl;            ///< the output level
	switch (outLvl) {
		case TrNone:
		case TrFatal:
		case TrError:
		case TrWarn:
		case TrInfo:
			hkwOutLvl = 0; // no output
			break;
		case TrInfo2:
			hkwOutLvl = 1; // + convergence info per trial
			break;
		case TrInfo3:
			hkwOutLvl = 2; // + convergence info per iter
			break;
		case TrDetail:
			hkwOutLvl = 3; // + info about cubic transf
			break;
		case TrDetail2:
			hkwOutLvl = 6; // + target matrices + decomp (hkwOutLvl=4)
			               // + Cholesky check
			break;
		case TrDetail3:
			hkwOutLvl = 7; // + cubic transf. detailed info
			break;
		case TrAll:
			hkwOutLvl = 11; // + matrix of outcomes
			break;
		default:
			std::cerr << "Warning: unknown value of 'outLvl' in "
			          << "ScenGen1PerHKW::gen_scen_main()" << std::endl;
			hkwOutLvl = 2;
			break;
	}
	// the output indicators
	double hkwOutErrMom;
	double hkwOutErrCorr;
	int hkwOutTrial;
	int hkwOutNmbIters;
	int hkwRetCode = ExtScenGenHKW::scengen_HKW (
		hkwTgMoms,  // target moments
		formOfMoms, // format of moments
		hkwTgCorrs, // target correlations
		prob,       // target probabilities (can be NULL)
		1,          // number of random variables
		nSc,        // number of scenarios to generate
		hkwOutSc,   // output - the array of generated scenario values
		0.01,       // maximal allowed error of moments
		0.01,       // maximal allowed error of correlations
		hkwOutLvl,  // level of output; computed from outLvl, see above
		2,          // max number of trials
		10,         // max number of iterations of the moment-matching heur.
		0,          // 0/1: use current values in hkwOutSc as start values?
		&hkwOutErrMom,  // output: final error of moments
		&hkwOutErrCorr, // output: final error of correlations
		&hkwOutTrial,   // output: number of used trials
		&hkwOutNmbIters // output: number of iterations (in the last trial)
	);
	if (hkwRetCode != 0) {
		MSG (TrError, "WARNING: scengen_HKW() finished with return code "
		              << hkwRetCode << " and total error of "
		              << hkwOutErrMom + hkwOutErrCorr << "!");
	} else {
		MSG (TrDetail, "scengen_HKW() finished with return code "
		               << hkwRetCode << " and total error of "
		               << hkwOutErrMom + hkwOutErrCorr);
	}

	// now sort the array
	std::sort(sortedVals.begin(), sortedVals.end());
}
