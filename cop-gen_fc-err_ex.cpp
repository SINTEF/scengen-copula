// examples of using the forecast-error-based generator
// NB: NOT USED ANYWHERE
//  - move some parts to test and delete the file!

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

struct BrInfo {
	int fromSc;
	int brTime;
};

class SC {
private:
	VectorI nOffspring;
	int nextSc;

	void count_offspring () {
		int T = branching.size();
		nOffspring.resize(T + 1);
		nOffspring(T) = 0;
		nOffspring(T - 1) = branching(T - 1);
		for (int t = T - 2; t >= 0; --t)
			nOffspring(t) = nOffspring(t + 1) * branching(t);
	}

public:
	VectorI branching;

	std::vector<BrInfo> get_struct_smps() {
		count_offspring();
		int nSc = nOffspring(0);
		auto brInfo = std::vector<BrInfo>(nSc);

		nextSc = 0;
		smps(0, 0, 0, 0, brInfo);

		return brInfo;
	}

	void smps(int curT, int curSc, int brT, int fromSc, std::vector<BrInfo> & brInfo) {
		if (nOffspring(curT) <= 1) {
			// no branching after this node
			brInfo[curSc].fromSc = fromSc;
			brInfo[curSc].brTime = brT;
		} else {
			// node has more than one offspring
			// first child - the same scenario
			smps(curT + 1, curSc, brT, fromSc, brInfo);
			// now other child (if any) - these are new scenarios
			for (int ch = 1; ch < branching(curT); ++ch) {
				nextSc++;
				smps(curT + 1, nextSc, curT + 1, curSc, brInfo);
			}
		}
	}

	std::vector<std::vector<BrInfo>> get_struct_txy() {
		count_offspring();
		int nSc = nOffspring(0);
		auto brInfo = std::vector<std::vector<BrInfo>>(nSc);

		nextSc = 0;
		int T = branching.size();
		VectorI curPath = ublas::zero_vector<int>(T+1);
		std::vector<VectorI> pathSc(nSc);
		txy_struct(0, curPath, pathSc);
		BrInfo branch;

		for (int sc = 0; sc < nSc; ++sc) {
			std::cout << "scenario " << sc << " : " << pathSc[sc] << '\n';
			branch.brTime = 0;
			branch.fromSc = pathSc[sc](0);
			brInfo[sc].push_back(branch);
			std::cout << brInfo[sc].size() << '\n';
			std::cout << brInfo[sc][0].brTime << "," << brInfo[sc][0].fromSc << '\n';
			for (int t = 1; t <= T; ++t) {
				if (pathSc[sc](t) != pathSc[sc](t-1)) {
					branch.brTime = t;
					branch.fromSc = pathSc[sc](t);
					brInfo[sc].push_back(branch);
				}
			}
		}
		return brInfo;
	}

	void txy_struct(int curT, VectorI & curPath, std::vector<VectorI> & pathSc) {
		if (nOffspring(curT) <= 1) {
			// no branching after this node
			int curSc = curPath(curPath.size() - 1);
			pathSc[curSc] = curPath;
		} else {
			// node had more than one offspring
			// first child - the same scenario
			txy_struct(curT + 1, curPath, pathSc);
			// now other child (if any) - these are new scenarios
			for (int ch = 1; ch < branching(curT); ++ch) {
				nextSc++;
				for (unsigned t = curT + 1; t < curPath.size(); ++t)
					curPath(t) = nextSc;
				txy_struct(curT + 1, curPath, pathSc);
			}
		}
	}
};

int main(int argc, char *argv[]) {

	SC scenTree;

	scenTree.branching.resize(8);
	scenTree.branching <<= 1, 1, 1, 3, 1, 1, 2, 1;

	auto treeStruct = scenTree.get_struct_smps();
	for (unsigned sc = 0; sc < treeStruct.size(); ++sc)
		std::cout << "scenario " << sc + 1 << " branches from "
		          << treeStruct[sc].fromSc + 1 << " at time "
		          << treeStruct[sc].brTime << "." << '\n';

	std::vector<std::vector<BrInfo>> treeStructTxy = scenTree.get_struct_txy();
	for (unsigned sc = 0; sc < treeStruct.size(); ++sc) {
		std::cout << "scenario " << sc + 1 << '\n';
		for (unsigned br = 0; br < treeStructTxy[sc].size(); ++br) {
			std::cout << " : " << treeStructTxy[sc][br].brTime << " -> "
			          <<treeStructTxy[sc][br].fromSc + 1 << '\n';
		}
	}

/*
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
		errors_to_scens(errScens, forecast, scens);
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
		                                false, perVarDt, intVarDt);
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
		                                false, perVarDt, intVarDt);
		FcErr_Gen::ScenTree scTree;

		VectorI branching(T);
		assert (T == 4 && "this test is fixed to four periods");
		branching <<= 3, 3, 2, 1;

		scenGen.gen_reg_tree(forecast, branching, scTree);
		scTree.display_per_scen();
		std::cout << std::endl;
		scTree.display_per_var();
	}
*/

	return 0;
}
