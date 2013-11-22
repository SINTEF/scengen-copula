/// \file cop-gen_lib.cpp
/**
	This file implements the interface of the cop-gen library
**/

#include "cop-gen_lib.hpp"

#include "copula-sample.hpp"
#include "margins.hpp"

using namespace ScenGenCop;


// output level - has to be defined
#ifdef NDEBUG
	OutputLevel defOutLvl = TrInfo;    // default value for release code
#else
	OutputLevel defOutLvl = TrDetail2; // default value for debug code
#endif
OutputLevel outLvl = defOutLvl;       // so it can be changed later


int gen_scen_normal(ublas::vector<double> const & tgMean,
                    ublas::vector<double> const & tgStD,
                    ublas::matrix<double> const & tgCorr,
                    unsigned const nSc,
                    ublas::matrix<double> & scens)
{
	auto p2tgCop = boost::make_shared<CopulaDef::CopInfoNormal>(tgCorr);
	auto p2tgMargins
		= boost::make_shared<MarginDistrib::NormalMargins>(tgMean, tgStD);

	// generate scenarios for the copula
	// TO DO: add possibility for shuffling?
	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	copSc.gen_sample();

	// transform margins to the target distribution
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks); // get the results in terms of ranks
	p2tgMargins->get_margin_distr(copRanks, scens); // convert to the target

	return 0;
}
