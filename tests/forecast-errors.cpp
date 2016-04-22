//#include "../cop-gen_lib.hpp"
#include "../copula-sample.hpp"
#include "../margins.hpp"
#include "../cop-gen_fc-err.hpp"

//#include <boost/numeric/ublas/assignment.hpp>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

using namespace Copula2D;

TEST(ForecastErrorsGen, InputData) {
	MatrixD histData;

	DimT N = 2;   // number of variables
	DimT T = 4;   // number of periods to generate
	DimT S = 12;  // number of scenarios to generate
	DimT D = 100; // number of data points in the file
	DimT i, r, s, t;

	// read the input file
	std::string histDataFName = "hist-fcasts.dat";
	std::ifstream histDataStr(histDataFName);
	if (!histDataStr) {
		throw std::ios_base::failure("Could not open input file `"
		                             + histDataFName + "'!");
	}
	histDataStr >> histData;
	//DISPLAY_NL(histData);
	EXPECT_EQ(histData.size1(), D) << "checking number of rows";
	EXPECT_EQ(histData.size2(), N * (T+1)) << "checking number of columns";

	// create a forecast equal to the first T observations
	MatrixD forecast(T, N);
	for (t = 0; t < T; ++t)
		for (i = 0; i < N; ++i)
			forecast(t, i) = histData(t, (T+1) * i);
	//DISPLAY_NL(forecast);

	// parameters for the generation of 2D copulas
	int perVarDt = 1;
	int intVarDt = 0;

	auto p2tgCop
		= boost::make_shared<FcErr_Gen::CopInfoForecastErrors>(
			N, histData, perVarDt, intVarDt);
	EXPECT_EQ(p2tgCop->dim(), N * T) << "checking dimension of the target cop.";
	EXPECT_EQ(p2tgCop->get_nmb_2d_copulas(), N*(T-1) + N*(N-1)/2*T)
		<< "checking the number of 2D targets"; // for perVarDt=1, intVarDt=0
	EXPECT_EQ(p2tgCop->hist_forecast_errors().size1(), N * T)
		<< "number of variables in the computed error matrix";
	EXPECT_EQ(p2tgCop->hist_forecast_errors().size2(), D - T)
		<< "number of data points in the computed error matrix";
	//DISPLAY_NL(p2tgCop->hist_forecast_errors());

	// NB: the margins-object needs the matrix of errors (used for generation),
	//     not the original matrix of historical data!
	auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>
		(p2tgCop->hist_forecast_errors());
	CopulaScen::CopulaSample copSc(p2tgCop, S);
	copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);
	EXPECT_EQ(copRanks.size1(), N * T) << "checking the number of variables";
	EXPECT_EQ(copRanks.size2(), S) << "checking the number of scenarios";

	// transform margins to the target distribution
	MatrixD errScens;  // dim = nVars, nSc]
	p2tgMargs->get_margin_distr(copRanks, errScens);
	EXPECT_EQ(errScens.size1(), N * T) << "checking the number of variables";
	EXPECT_EQ(errScens.size2(), S) << "checking the number of scenarios";
	//DISPLAY_NL(errScens);

	// convert to the scenarios for the original multi-period problem
	std::vector<MatrixD> scens;  // dim = nSc * [T, N]
	FcErr_Gen::errors_to_values(errScens, forecast, scens);
	//for (s = 0; s < S; ++s) {
	//	DISPLAY_NL(scens[s]);
	//}

	for (s = 0; s < S; ++s) {
		for (t = 0; t < T; ++t) {
			for (i = 0; i < N; ++i) {
				r = N * t + i; // rows grouped by time
				EXPECT_DOUBLE_EQ(scens[s](t, i),
				                 forecast(t, i) + errScens(r, s))
				                 << "checking value (" << s << ", " << i
				                 << ", " << t << ")";
			}
		}
	}
}
