/// \file cop-gen_lib.cpp
/**
	This file implements the interface of the cop-gen library
**/

#include "cop-gen_lib.hpp"
#include "copula-sample.hpp"
#include "margins.hpp"

// output/verbosity level - must exist; the value can be changed
OutputLevel outLvl = static_cast<OutputLevel>(ScenGenCop::defVerb);


// generate scenarios from normal distribution
double ScenGenCop::gen_scen_normal(ublas::vector<double> const & tgMean,
                       ublas::vector<double> const & tgStD,
                       ublas::matrix<double> const & tgCorr,
                       unsigned const nSc,
                       ublas::matrix<double> & scens,
                       bool const varInC,
                       unsigned const verb)
{
	if (verb > OutputLevel::TrAll)
		outLvl = OutputLevel::TrAll;
	else
		outLvl = static_cast<OutputLevel>(verb);

	auto p2tgCop = boost::make_shared<CopulaDef::CopInfoNormal>(tgCorr);
	auto p2tgMargins
		= boost::make_shared<MarginDistrib::NormalMargins>(tgMean, tgStD);

	// generate scenarios for the copula
	// TO DO: add possibility for shuffling?
	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	double scCopDist = copSc.gen_sample();

	// transform margins to the target distribution
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks); // get the results in terms of ranks

	// convert to the target distribution
	if (varInC) {
		// want variables in columns -> have to transpose the results
		MatrixD trScens;
		p2tgMargins->get_margin_distr(copRanks, trScens);
		scens = ublas::trans(trScens);
	} else {
		p2tgMargins->get_margin_distr(copRanks, scens);
	}

	return scCopDist;
}


// generate scenarios based on a sample from the distrib. (historical data)
double ScenGenCop::gen_scen_sample(
	ublas::matrix<double> const & tgSample,
	unsigned const nSc, ublas::matrix<double> & scens,
	bool const inMatPerC, bool const outMatPerC,
	unsigned const verb)
{
	if (verb > OutputLevel::TrAll)
		outLvl = OutputLevel::TrAll;
	else
		outLvl = static_cast<OutputLevel>(verb);

	// the constructors below expect a row-wise target matrix ([nVar * nPts])
	ublas::matrix<double> tgSampT;
	if (inMatPerC)
		tgSampT = ublas::trans(tgSample);
	ublas::matrix<double> const & ref2tg = (inMatPerC ? tgSampT : tgSample);
	auto p2tgCop = boost::make_shared<CopulaDef::CopInfoData>(ref2tg);
	auto p2tgMargins = boost::make_shared<MarginDistrib::SampleMargins>(ref2tg);

	// generate scenarios for the copula
	// TO DO: add possibility for shuffling?
	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	double scCopDist = copSc.gen_sample();

	// transform margins to the target distribution
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks); // get the results in terms of ranks

	// convert to the target distribution
	if (outMatPerC) {
		// want variables in columns -> have to transpose the results
		MatrixD trScens;
		p2tgMargins->get_margin_distr(copRanks, trScens);
		scens = ublas::trans(trScens);
	} else {
		p2tgMargins->get_margin_distr(copRanks, scens);
	}

	return scCopDist;
}
