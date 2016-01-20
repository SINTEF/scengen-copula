#include "cop-gen_fc-err.hpp"

#include <iostream>
#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>
//#include <cassert>

//using namespace std;
using namespace CopulaDef;
using namespace Copula2D;
using namespace FcErr_Gen;


// --------------------------------------------------------------------------
// static non-class functions - just for internal use

/// get row of forecast for var. \c i and step \c dt in the historical data
static DimT hist_data_row_of(DimT const i, DimT const dt,
                             DimT const N, DimT const T)
	{ return (dt - 1) * N + i; }

static auto firstRows(MatrixD const & X, DimT nR) -> decltype(ublas::subrange(X, 0, 0, 0, 0))
{
	return ublas::subrange(X, 0, nR, 0, X.size2());
}


// --------------------------------------------------------------------------
// non-class functions from the FcErr_Gen namespace

/*
	\param[in]  histData  historical values and forecasts [nPts, N * (T+1)];
	                      columns are grouped by variables, for each (i, t)
	                      we have: val(i,t), fc(i,t,t+1), fc(i,t,t+2), ...
	                      where fc(i,t,t+dt) is forecast made at t for t+dt
	\param[out] histFErr  the computed historical forecast errors [nVars, nPts]
	\param[in]         N  the original dimension (number of orig. variables)
	\param[in]  dataSort  sorting/structure of \c histData
*/
void FcErr_Gen::histdata_to_errors(MatrixD const & histData, MatrixD & histFErr,
                                   DimT const N, HistDataSort const dataSort)
{
	DimT T = histData.size2() / N - 1;
	if (histData.size2() != N * (T+1))
		throw std::length_error("inconsistence in the size of historical data");
	DimT nPts = histData.size1() - T; // need extra data to compute the errors
	DimT nVars = N * T;

	DimT i, dt, cF, cV, rE, j;

	// fill the histFErr matrix (NB: transposed rel. to histData!)
	histFErr.resize(nVars, nPts);
	for (i = 0; i < N; ++i) {
		cV = i * (T+1); // column of the variable value
		for (dt = 1; dt <= T; ++dt) {
			cF = cV + dt;      // column of the forecast
			rE = hist_data_row_of(i, dt, N, T); // row within histFErr
			for (j = 0; j < nPts; ++j) {
				switch (dataSort) {
				case HistDataSort::fCastTimeAsc:
					histFErr(rE, j) = histData(j, cF) - histData(j + dt, cV);
					break;
				default:
					throw std::logic_error("unsupported sorting option");
				}
			}
		}
	}
}


// convert scenarios of errors to scenarios of the original values
/*
	\param[in]  errSc  error-scenarios; [nVars, nSc]
	\param[in] forecast  forecast for the whole time horizon; [T, N]
	\param[out] scens  output scenarios; nSc * [T, N]
*/
void FcErr_Gen::errors_to_values(MatrixD const & errSc, MatrixD const & forecast,
                                 std::vector<MatrixD> & scens)
{
	DimT nVars = errSc.size1();
	DimT nSc = errSc.size2();
	DimT N = forecast.size2();
	DimT T = nVars / N;
	if (nVars != T * N)
		throw std::length_error
			("inconsistent sizes of error and forecast matrices");
	if (forecast.size1() < T)
		throw std::length_error
			("not enough periods in the forecast matrix");

	DimT i, r, s, t;

	scens.resize(nSc);
	for (s = 0; s < nSc; ++s) {
		scens[s].resize(T, N);
		for (t = 0; t < T; ++t) {
			for (i = 0; i < N; ++i) {
				r = N * t + i; // row in errSc
				scens[s](t, i) = forecast(t, i) + errSc(r, s);
			}
		}
	}
}


// --------------------------------------------------------------------------
// class CopInfoForecastErrors

// constructor with the historical forecast errors
/*
	\param[in]  nmbVars  number of variables
	\param[in] histFErr  historical forecast errors [N * T, nPts];
	                     columns are grouped by periods: val(i,dt) for
	                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
	\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (i, dt ± u) for u = 1..perVarDt
	\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (j, dt ± u) for u = 0..intVarDt
*/
CopInfoForecastErrors::CopInfoForecastErrors(DimT const nmbVars,
	                                         MatrixD const & histFErr,
	                                         int perVarDt, int intVarDt)
: CopInfoData(histFErr),
  N(nmbVars), T(histFErr.size1() / N)
{
	if (histFErr.size1() != N * T)
		throw std::length_error("inconsistent dimensions of input data");

	setup_2d_targets(perVarDt, intVarDt);
}


// constructor with the historical forecast errors
/*
	\param[in]  nmbVars  number of variables
	\param[in] histData  historical values and forecasts [nPts, N * (T+1)];
	                     columns are grouped by variables, for each (i, t)
	                     we have: val(i,t), fc(i,t,t+1), fc(i,t,t+2), ...
	                     where fc(i,t,t+dt) is forecast made at t for t+dt
	\param[in] dataSort  sorting/structure of \c histData
	\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (i, dt ± u) for u = 1..perVarDt
	\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (j, dt ± u) for u = 0..intVarDt
*/
CopInfoForecastErrors::CopInfoForecastErrors(DimT const nmbVars,
                                             MatrixD const & histData,
                                             HistDataSort const dataSort,
                                             int perVarDt, int intVarDt)
: CopInfoData(histData.size2() - nmbVars, true), // the first param is N * T
  N(nmbVars), T(histData.size2() / nmbVars - 1)
{
	if (histData.size2() != N * (T+1))
		throw std::length_error("inconsistent dimensions of input data");

	// set the number of historical data points (defined in CopInfoData)
	// - smaller than the the size of histData, because we lose some
	//   data for computing the errors/differences
	nPts = histData.size1() - T;

	// fill the hData matrix (NB: transposed rel. to histData!)
	histdata_to_errors(histData, hData, N, dataSort);

	// rest of the setup
	fill_ranks_etc();                     // fill hRanks and hU01
	setup_2d_targets(perVarDt, intVarDt); // create the bivariate copula objects
}


// version of \c hist_data_row_of() for use inside the class (fewer parameters)
DimT CopInfoForecastErrors::row_of(DimT const i, DimT const dt)
{
	return hist_data_row_of(i, dt, N, T);
}



// creates objects for the 2D targets; called from the constructors
/*
	\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (i, dt ± u) for u = 1..perVarDt
	\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (j, dt ± u) for u = 0..intVarDt
*/
void CopInfoForecastErrors::setup_2d_targets(int perVarDt, int intVarDt)
{
	assert (nVars > 0 && "must have a known dimension by now");
	if (p2Info2D.size1() != nVars || p2Info2D.size2() != nVars) {
		if (p2Info2D.size1() * p2Info2D.size2() > 0) {
			std::cerr << "Warning: resizing the matrix of 2D-copula objects!"
			          << std::endl;
		}
		p2Info2D.resize(nVars, nVars);
	}

	DimT i, j, v1, v2;
	int dt, u;

	for (i = 0; i < N; ++i) {
		for (dt = 1; dt <= (int) T; ++dt) {
			v1 = row_of(i, dt);
			assert (v1 >= 0 && v1 < nVars && "bound check");
			// add copulas for one variable and varying dt
			for (u = 1; u <= perVarDt; ++u) {
				if (dt - u >= 1) {
					v2 = row_of(i, dt - u);
					assert (v2 >= 0 && v1 < nVars && "bound check");
					if (v2 > v1)
						attach_2d_target(new Cop2DData(hData, v1, v2, this),
						                 v1, v2);
				}
				if (dt + u <= (int) T) {
					v2 = row_of(i, dt + u);
					assert (v2 >= 0 && v1 < nVars && "bound check");
					if (v2 > v1)
						attach_2d_target(new Cop2DData(hData, v1, v2, this),
						                 v1, v2);
				}
			}
			// add copulas between variables and possibly varying dt
			for (j = 0; j < N; ++j) {
				if (j != i) {
					for (u = 0; u <= intVarDt; ++u) {
						if (dt - u >= 1) {
							v2 = row_of(j, dt - u);
							assert (v2 >= 0 && v1 < nVars && "bound check");
							if (v2 > v1)
								attach_2d_target(
									new Cop2DData(hData, v1, v2, this), v1, v2);
						}
						if (dt + u <= (int) T) {
							v2 = row_of(j, dt + u);
							assert (v2 >= 0 && v1 < nVars && "bound check");
							if (v2 > v1)
								attach_2d_target(
									new Cop2DData(hData, v1, v2, this), v1, v2);
						}
					}
				}
			}
		}
	}
}


// --------------------------------------------------------------------------
// class ScenTree

bool ScenTree::is_empty() const
{
	return (scenVals.size() == 0);
}

/// get the number of periods
DimT ScenTree::nmb_periods() const
{
	if (scenVals.size() == 0)
		return 0;
	return scenVals[0].size1(); // assuming all scenarios are equal
}

/// get the number of variables
DimT ScenTree::nmb_variables() const
{
	if (scenVals.size() == 0)
		return 0;
	return scenVals[0].size2(); // assuming all scenarios are equal
}


void ScenTree::display_per_scen() const
{
	DimT s, S = scenVals.size();
	for (s = 0; s < S; ++s) {
		std::cout << "values for scenario " << s+1 << ":" << std::endl
		          << scenVals[s] << std::endl;
	}
}

void ScenTree::display_per_var() const
{
	DimT s, S = scenVals.size();
	DimT t, T = scenVals[0].size1();
	DimT i, N = scenVals[0].size2();
	for (i = 0; i < N; ++i) {
		std::cout << "values for margin/variable " << i+1 << ":" << std::endl;
		if (rootVals.size() == N) {
			// copy the root value to all scenarios
			for (s = 0; s < S; ++s)
				std::cout << rootVals[i] << "\t";
			std::cout << std::endl;
		}

		for (t = 0; t < T; ++t) {
			for (s = 0; s < S; ++s)
				std::cout << scenVals[s](t, i) << "\t";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}


// --------------------------------------------------------------------------
// class FcErrTreeGen


// constructor with the historical values and forecasts
/*
	\param[in]  nmbVars  number of stoch. variables
	\param[in] histData  historical values and forecasts [nPts, N * (T+1)];
	                     columns are grouped by variables, for each (i, t)
	                     we have: val(i,t), fc(i,t,t+1), fc(i,t,t+2), ...
	                     where fc(i,t,t+dt) is forecast made at t for t+dt
	\param[in] dataSort  sorting/structure of \c histData
	\param[in] maxPerVarDt  for var. (i, dt), we generate 2D-copulas with
	                        (i, dt ± u) for u = 1..perVarDt
	\param[in] maxIntVarDt  for var. (i, dt), we generate 2D-copulas with
	                        (j, dt ± u) for u = 0..intVarDt
*/
FcErrTreeGen::FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
                           HistDataSort const dataSort,
                           int maxPerVarDt, int maxIntVarDt)
: N(nmbVars), T(histData.size2() / N - 1),
  perVarDt(maxPerVarDt), intVarDt(maxIntVarDt)
{
	if (histData.size2() != N * (T+1))
		throw std::length_error("inconsistent dimensions of input data");

	// create the matrix of historical forecast errors
	// note: We cannot use histFErr here, since it is a ref. to const matrix.
	//       However, it works with _histFErr, which is a non-const matrix
	//       .. even if histFErr references it .. weird but true
	histdata_to_errors(histData, _histFErr, N, dataSort);
}

// generate a 2-stage tree (a fan), with given number of scenarios
/*
	\param[in] forecast  forecast for the whole horizon; [T, N]
	\param[in]    nScen  number of scenarios to generate
	\param[out] outTree  the resulting scenario tree
*/
void FcErrTreeGen::gen_2stage_tree(MatrixD const & forecast, DimT const nScen,
                                   ScenTree & outTree)
{
	auto p2tgCop = boost::make_shared<CopInfoForecastErrors>(
		N, histFErr, perVarDt, intVarDt);

	auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>
		(histFErr);

	CopulaScen::CopulaSample copSc(p2tgCop, nScen);
	copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);

	// transform margins to the target distribution
	MatrixD errScens;  // dim = nVars, nSc]
	p2tgMargs->get_margin_distr(copRanks, errScens);
	//DISPLAY_NL(errScens);

	// convert to the scenarios for the original multi-period problem
	errors_to_values(errScens, forecast, outTree.scenVals);
}


// add one stage to an existing scenario tree
/*
	\param[in]         nBr  number of branches to create at each scenario
	\param[in]          dT  number of periods to add to the tree (>=1)
	\param[out] p2copRanks  the generated copula ranks, [nVar, nSc]
	\param[out]  totErrors  converted output values, [T * N, nSc]
	\param[in] p2prevRanks  copula ranks from the prev. iteration, [nVar, nSc]
	\param[in]  p2forecast  vector of forecasts for the whole horizon; [T, N]
	\param[out]  p2outTree  the resulting scenario tree (if required)

	// the output scenario tree is needed only at the end
	// in previous iterations, the last parameters would be empty
*/
void FcErrTreeGen::add_one_stage(DimT const nBr, DimT const dT,
                                 MatrixI const & prevRanks,
                                 MatrixI & copRanks, MatrixD & totErrors)
{
	DimT nPers = dT; // number of periods to generate
	DimT nSc = nBr;  // number of scenarios to generate
	DimT nFixM = 0;  // number of margins to fix in the new copula
	DimT i, s;

	// we use the fact that the columns of histFErr are sorted by periods:
	// first N variables of the 1st period, then the 2nd period, etc.

	if (prevRanks.size1() > 0) {
		// have a non-empty matrix of existing copula (from previous iteration)
		nFixM = prevRanks.size1();
		nPers += nFixM / N; // have N margins per period
		nSc *= prevRanks.size2();
	}

	MatrixD const & errM = (nPers == T ? histFErr
	                                   : firstRows(histFErr, nPers * N));
	MSG(TrDetail2, "historical data for the new stage:" << std::endl << errM);

	auto p2tgCop = boost::make_shared<CopInfoForecastErrors>(
		N, errM, perVarDt, intVarDt);

	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	// fixed the starting copula, if we have one (expanded to nSc)
	if (nFixM > 0) {
		MatrixI fixedMargs(nFixM, nSc);
		DimT s0, br;
		s = 0;
		for (s0 = 0; s0 < prevRanks.size2(); ++s0) {
			for (br = 0; br < nBr; ++br) {
				// cannot work with whole columns, since boost does not have
				// vector + scalar overloaded
				for (i = 0; i < nFixM; ++i)
					fixedMargs(i, s) = nBr * prevRanks(i, s0) + nBr / 2;
				//ublas::column(fixedMargs, s) = ublas::column(*p2prevRanks, s0);
				++s;
			}
		}
		assert (s == nSc && "consistency check");
		DBGSHOW_NL(TrDetail, prevRanks);
		DBGSHOW_NL(TrDetail, fixedMargs);
		copSc.fix_marg_values(fixedMargs, true);
	}

	copSc.gen_sample();
	copSc.get_result_ranks(copRanks);
	DBGSHOW_NL(TrDetail, copRanks);

	// transform margins to the target distribution
	MatrixD errScens;  // dim = nVars, nSc]
	auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>(errM);
	p2tgMargs->get_margin_distr(copRanks, errScens);
	MSG(TrDetail, "generated matrix of errors:" << std::endl << errScens);

	// now, copy the result to totErrors, but only the periods generated here
	if (totErrors.size1() > 0) {
		DimT scenMult = totErrors.size2() / nSc;
		if (nSc * scenMult != totErrors.size2())
			throw std::length_error("inconsistent size of the output matrix");
		for (i = nFixM; i < nFixM + N * dT; ++i) {
			for (s = 0; s < nSc; ++s)
				for (DimT sM = 0; sM < scenMult; ++sM)
					totErrors(i, scenMult * s + sM) = errScens(i, s);
		}
	}

	/*
	// check if the caller has asked to create the scenario tree
	if (p2outTree) {
		if (!p2forecast)
			throw std::logic_error("non-zero outTree, but empty forecast");

		auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>(errM);

		// transform margins to the target distribution
		MatrixD errScens;  // dim = [nVars, nSc]
		p2tgMargs->get_margin_distr(*p2copRanks, errScens);
		DBGSHOW(TrDetail, errScens);

		// convert to the scenarios for the original multi-period problem
		// note: we cannot write p2outTree->branching, as we do not know
		//       about branching in the previous periods .. left to the caller
		errors_to_values(errScens, *p2forecast, p2outTree->scenVals);
	}
	*/
}


// generate a regular multistage tree, with given branching per period
/*
	\param[in]  forecast  forecast for the whole horizon; [T, N]
	\param[in] branching  branching factor per period
	\param[out]  outTree  the resulting scenario tree

	\note By a regular tree we mean a scenario tree where all nodes
	      in the same period have the same number of child nodes.
*/
void FcErrTreeGen::gen_reg_tree(MatrixD const & forecast,
                                VectorI const & branching, ScenTree & outTree)
{
	if (forecast.size2() != N)
		throw std::length_error("wrong number of variables in the forecast");
	DimT nmbPer = T;  // number of periods we generate scenarios for
	if (branching.size() != forecast.size1())
		WARNING("different number of periods in branching vector and forecast")
	if (forecast.size1() < nmbPer)
		nmbPer = forecast.size1();
	if (branching.size() < nmbPer)
		nmbPer = branching.size();

	MatrixI copRanks[2];
	DimT outCopIndx = 0; // index of the output copula out of the two
	//MatrixI copRanks1, copRanks2;
	//MatrixI & initCopR = copRanks2;
	//MatrixI & genCopR = copRanks1;

	DimT t, nSc = 1;  // total number of scenarios
	for (t = 0; t < nmbPer; ++t)
		nSc *= branching(t);
	MatrixD errScens  // overall output matrix, initialized to zeros
		= ublas::zero_matrix<double>(N * nmbPer, nSc);
	DimT tBr = 0;     // period of the current branching
	DimT tBrNext;     // period of the next branching
	do {
		INFO("generating subtree starting at period " << tBr + 1);

		// find the next branching period
		tBrNext = tBr + 1;
		while (tBrNext < nmbPer && branching(tBrNext) == 1)
			++tBrNext;

		// not the last branching
		add_one_stage(branching(tBr), tBrNext - tBr,
		              copRanks[1-outCopIndx], copRanks[outCopIndx], errScens);
		MSG(TrDetail, "updated overall matrix of errors:" << std::endl
		              << errScens);

		outCopIndx = 1 - outCopIndx; // swap input and output copulas
		tBr = tBrNext;
	} while (tBrNext < nmbPer);

	DBGSHOW_NL(TrDetail, forecast);
	errors_to_values(errScens, forecast, outTree.scenVals);
	outTree.branching = branching;
}


// --------------------------------------------------------------------------
// methods

// add one stage to an existing scenario tree
/*
	\param[in] initTree  the tree we are adding to (might be empty)
	\param[out] outTree  the extended tree
	\param[in]    tgCop  specifications for generating the target cop
	\param[in]      nBr  number of branches to create at each scenario
	\param[in]       dT  number of periods to add to the tree (>=1)
*//*
void add_one_stage(ScenTree const & initTree, ScenTree & outTree,
                   CopulaDef::CopInfoForecastErrors const & tgCop,
                   DimT const nBr, DimT const dT = 1)
{
	if (initTree.is_empty()) {

	}
}
*/
