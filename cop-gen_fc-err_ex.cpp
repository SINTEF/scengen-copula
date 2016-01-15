/// main function for the forecast-error-based generator
/**
	\note The main function cannot be in \c cop-gen_fc-err.cpp, since
	      this file is used also with other build targets.
**/

#include "cop-gen_fc-err.hpp"
#include "copula-sample.hpp"
#include "margins.hpp"

#include <boost/numeric/ublas/assignment.hpp>
#include <iostream>
#include <fstream>

using namespace Copula2D;
using namespace FcErr_Gen;
//using DataSort = CopulaDef::CopInfoForecastErrors::HistDataSort;

// output level
#ifdef NDEBUG
	OutputLevel defOutLvl = TrInfo;    // default value for release code
#else
	OutputLevel defOutLvl = TrDetail2; // default value for debug code
#endif
OutputLevel outLvl = defOutLvl; // can be changed later


int main(int argc, char *argv[]) {
	ECHO ("Example 1: simple 2-stage tree, built manually");
	{
		DimT N; // number of variables
		DimT T; // number of periods to generate
		DimT S; // number of scenarios to generate
		DimT i, s, t;

		MatrixD histData;
		// read the input file
		std::string histDataFName = "test_fcast-err.dat";
		std::ifstream histDataStr(histDataFName);
		if (!histDataStr) {
			throw std::ios_base::failure("Could not open input file `"
			                             + histDataFName + "'!");
		}
		histDataStr >> histData;
		histDataStr.close();

		// Normally, N and T would be inferred from the forecast matrix
		//  - the only thing we get here is: histData.size2() = N * (T+1)
		// Here, we have to fake forecast, so we have to fixed them here:
		N = 2;
		T = histData.size2() / N - 1;
		if (N * (T+1) != histData.size2())
			throw std::domain_error("wrong input dimensions");
		//
		// create a forecast equal to the first T observations
		MatrixD forecast(T, N);
		for (t = 0; t < T; ++t)
			for (i = 0; i < N; ++i)
				forecast(t, i) = histData(t, (T+1) * i);
		//DISPLAY_NL(forecast);

		// parameters for the generation of 2D copulas
		int perVarDt = 1;
		int intVarDt = 0;
		S = 12;

		auto p2tgCop
			= boost::make_shared<CopInfoForecastErrors>(
				N, histData, HistDataSort::fCastTimeAsc, perVarDt, intVarDt);
		//DISPLAY_NL(p2tgCop->hist_forecast_errors());

		// NB: the margins-object needs the matrix of errors (used for generation),
		//     not the original matrix of historical data!
		auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>
			(p2tgCop->hist_forecast_errors());
		CopulaScen::CopulaSample copSc(p2tgCop, S);
		copSc.gen_sample();
		MatrixI copRanks;
		copSc.get_result_ranks(copRanks);

		// transform margins to the target distribution
		MatrixD errScens;  // dim = nVars, nSc]
		p2tgMargs->get_margin_distr(copRanks, errScens);
		//DISPLAY_NL(errScens);

		// convert to the scenarios for the original multi-period problem
		std::vector<MatrixD> scens;  // dim = nSc * [T, N]
		errors_to_values(errScens, forecast, scens);
		for (s = 0; s < S; ++s) {
			ECHO("values for scenario " << s+1 << ": " << std::endl << scens[s]);
		}
	}
	ECHO ("Example 1 finished." << std::endl);

	ECHO ("Example 2: the same 2-stage tree, using FcErrTreeGen");
	{
		DimT N; // number of variables
		DimT T; // number of periods to generate
		DimT S; // number of scenarios to generate
		DimT i, t;

		MatrixD histData;
		// read the input file
		std::string histDataFName = "test_fcast-err.dat";
		std::ifstream histDataStr(histDataFName);
		if (!histDataStr) {
			throw std::ios_base::failure("Could not open input file `"
			                             + histDataFName + "'!");
		}
		histDataStr >> histData;
		histDataStr.close();

		// Normally, N and T would be inferred from the forecast matrix
		//  - the only thing we get here is: histData.size2() = N * (T+1)
		// Here, we have to fake forecast, so we have to fixed them here:
		N = 2;
		T = histData.size2() / N - 1;
		if (N * (T+1) != histData.size2())
			throw std::domain_error("wrong input dimensions");
		//
		// create a forecast equal to the first T observations
		MatrixD forecast(T, N);
		for (t = 0; t < T; ++t)
			for (i = 0; i < N; ++i)
				forecast(t, i) = histData(t, (T+1) * i);
		//DISPLAY_NL(forecast);

		// parameters for the generation of 2D copulas
		int perVarDt = 1;
		int intVarDt = 0;
		S = 12;

		FcErr_Gen::FcErrTreeGen scenGen(N, histData, HistDataSort::fCastTimeAsc,
		                                perVarDt, intVarDt);
		FcErr_Gen::ScenTree scTree;
		scenGen.gen_2stage_tree(forecast, S, scTree);
		scTree.display_per_scen();
	}

	ECHO ("Example 3: a 3-stage tree, using FcErrTreeGen");
	{
		DimT N; // number of variables
		DimT T; // number of periods to generate
		//DimT S; // number of scenarios to generate
		DimT i, t;

		MatrixD histData;
		// read the input file
		std::string histDataFName = "test_fcast-err.dat";
		std::ifstream histDataStr(histDataFName);
		if (!histDataStr) {
			throw std::ios_base::failure("Could not open input file `"
			                             + histDataFName + "'!");
		}
		histDataStr >> histData;
		histDataStr.close();

		// Normally, N and T would be inferred from the forecast matrix
		//  - the only thing we get here is: histData.size2() = N * (T+1)
		// Here, we have to fake forecast, so we have to fixed them here:
		N = 2;
		T = histData.size2() / N - 1;
		if (N * (T+1) != histData.size2())
			throw std::domain_error("wrong input dimensions");
		//
		// create a forecast equal to the first T observations
		MatrixD forecast(T, N);
		for (t = 0; t < T; ++t)
			for (i = 0; i < N; ++i)
				forecast(t, i) = histData(t, (T+1) * i);
		//DISPLAY_NL(forecast);

		// parameters for the generation of 2D copulas
		int perVarDt = 1;
		int intVarDt = 0;
		//S = 12;

		FcErr_Gen::FcErrTreeGen scenGen(N, histData, HistDataSort::fCastTimeAsc,
		                                perVarDt, intVarDt);
		FcErr_Gen::ScenTree scTree;

		VectorI branching(T);
		assert (T == 4 && "this test is fixed to four periods");
		branching <<= 4, 1, 3, 1;

		scenGen.gen_reg_tree(forecast, branching, scTree);
		scTree.display_per_scen();
	}

	return 0;
}
