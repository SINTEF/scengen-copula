#include "cop-gen_fc-err.hpp"

#include <iostream>
#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>
//#include <cassert>

using namespace std;
using namespace CopulaDef;
using namespace Copula2D;
using namespace FcErr_Gen;

using DataSort = CopulaDef::CopInfoForecastErrors::HistDataSort;

// --------------------------------------------------------------------------
// class CopInfoForecastErrors

// constructor with the historical forecast errors
/*
	\param[in] forecast  forecast for the whole time horizon [N, T]
	\param[in] histFErr  historical forecast errors [N * T, nPts];
	                     columns are grouped by periods: val(i,dt) for
	                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
	\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (i, dt ± u) for u = 1..perVarDt
	\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
	                     (j, dt ± u) for u = 0..intVarDt
*/
CopInfoForecastErrors::CopInfoForecastErrors(MatrixD const & forecast,
	                                         MatrixD const & histFErr,
	                                         int perVarDt, int intVarDt)
: CopInfoData(histFErr),
  N(forecast.size1()), T(forecast.size2()), fCast(forecast)
{
	setup_2d_targets(perVarDt, intVarDt);
}


// constructor with the historical forecast errors
/*
	\param[in] forecast  forecast for the whole time horizon [N, T]
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
CopInfoForecastErrors::CopInfoForecastErrors(MatrixD const & forecast,
                                             MatrixD const & histData,
                                             HistDataSort const dataSort,
                                             int perVarDt, int intVarDt)
: CopInfoData(forecast.size1() * forecast.size2(), true),
  N(forecast.size1()), T(forecast.size2()), fCast(forecast)
{
	// set the number of historical data points (defined in CopInfoData)
	// - smaller than the the size of histData, because we lose some
	//   data for computing the errors/differences
	nPts = histData.size1() - T;

	// fill the hData matrix (NB: transposed rel. to histData!)
	histdata_to_errors(histData, hData, N, dataSort);
	/*
	hData.resize(nVars, nPts);
	for (i = 0; i < N; ++i) {
		cV = i * (T+1); // column of the variable value
		for (dt = 1; dt <= T; ++dt) {
			cF = cV + dt;      // column of the forecast
			rE = rowOf(i, dt); // row within hData
			for (j = 0; j < nPts; ++j) {
				switch (dataSort) {
				case HistDataSort::fCastTimeAsc:
					hData(rE, j) = histData(j, cF) - histData(j + dt, cV);
					break;
				default:
					throw std::logic_error("unsupported sorting option");
				}
			}
		}
	}
	*/

	// rest of the setup
	fill_ranks_etc();                     // fill hRanks and hU01
	setup_2d_targets(perVarDt, intVarDt); // create the bivariate copula objects
}


/*
	\param[in]  histData  historical values and forecasts [nPts, N * (T+1)];
	                      columns are grouped by variables, for each (i, t)
	                      we have: val(i,t), fc(i,t,t+1), fc(i,t,t+2), ...
	                      where fc(i,t,t+dt) is forecast made at t for t+dt
	\param[out] histFErr  the computed historical forecast errors [nVars, nPts]
	\param[in]         N  the original dimension (number of orig. variables)
	\param[in]  dataSort  sorting/structure of \c histData
*/
void CopInfoForecastErrors::histdata_to_errors(
	MatrixD const & histData, MatrixD & histFErr,
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
			rE = rowOf(i, dt, N, T); // row within histFErr
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
			cerr << "Warning: resizing the matrix of 2D-copula objects!" << endl;
		}
		p2Info2D.resize(nVars, nVars);
	}

	DimT i, j, v1, v2;
	int dt, u;

	for (i = 0; i < N; ++i) {
		for (dt = 1; dt <= (int) T; ++dt) {
			v1 = rowOf(i, dt);
			assert (v1 >= 0 && v1 < nVars && "bound check");
			// add copulas for one variable and varying dt
			for (u = 1; u <= perVarDt; ++u) {
				if (dt - u >= 1) {
					v2 = rowOf(i, dt - u);
					assert (v2 >= 0 && v1 < nVars && "bound check");
					if (v2 > v1)
						attach_2d_target(new Cop2DData(hData, v1, v2, this),
						                 v1, v2);
				}
				if (dt + u <= (int) T) {
					v2 = rowOf(i, dt + u);
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
							v2 = rowOf(j, dt - u);
							assert (v2 >= 0 && v1 < nVars && "bound check");
							if (v2 > v1)
								attach_2d_target(
									new Cop2DData(hData, v1, v2, this), v1, v2);
						}
						if (dt + u <= (int) T) {
							v2 = rowOf(j, dt + u);
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

// convert scenarios of errors to scenarios of the original values
/*
	\param[in]  errSc  error-scenarios; [nVars, nSc]
	\param[out] scens  output scenarios; nSc * [T, N]
*/
void CopInfoForecastErrors::errors_to_values(MatrixD const & errSc,
                                             std::vector<MatrixD> & scens) const
{
	if (errSc.size1() != nVars)
		throw std::length_error("wrong number of variables in errSc");
	DimT nSc = errSc.size2();

	DimT i, r, s, t;

	scens.resize(nSc);
	for (s = 0; s < nSc; ++s) {
		scens[s].resize(T, N);
		for (t = 0; t < T; ++t) {
			for (i = 0; i < N; ++i) {
				r = N * t + i; // row in errSc
				scens[s](t, i) = fCast(i, t) + errSc(r, s);
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

void ScenTree::display_per_scen() const
{
	DimT s, S = scenVals.size();
	for (s = 0; s < S; ++s) {
		ECHO("values for scenario " << s+1 << ": " << std::endl << scenVals[s]);
	}
}


// --------------------------------------------------------------------------
// class FcErrTreeGen


// constructor with the historical values and forecasts
/*
	\param[in] forecast  forecast for the whole time horizon [N, T]
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
FcErrTreeGen::FcErrTreeGen(MatrixD const & forecast, MatrixD const & histData,
                           DataSort const dataSort,
                           int maxPerVarDt, int maxIntVarDt)
: N(forecast.size1()), T(forecast.size2()),
  perVarDt(maxPerVarDt), intVarDt(maxIntVarDt), curFcast(forecast)
{
// create the matrix of historical data
	// note: since histFErr is a reference to _histFErr, we are actually
	//       writing the data to _histFErr
	CopulaDef::CopInfoForecastErrors::histdata_to_errors(
		histData, _histFErr, N, dataSort);
}

// generate a 2-stage tree (a fan), with given number of scenarios
/*
	\param[in] nScen  number of scenarios to generate
	\param[in]  nPer  number of periods to generate, nPer <= T; 0 means T
*/
void FcErrTreeGen::gen_2stage_tree(DimT const nScen, ScenTree & outTree,
                                   DimT const nPer)
{
	if (nPer != 0 && nPer != T) {
		// shorter horizon than what we have forecast for
		//  - not yet implemented!
		//  - would need either extra parameters to other methods,
		//    or, probably easier, make temp. submatrices with given length
		throw std::length_error("short horizons not yet implemented");
	}

	auto p2tgCop = boost::make_shared<CopulaDef::CopInfoForecastErrors>(
		curFcast, histFErr, perVarDt, intVarDt);

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
	//std::vector<MatrixD> scens;  // dim = nSc * [T, N]
	p2tgCop->errors_to_values(errScens, outTree.scenVals);
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
*/
void add_one_stage(ScenTree const & initTree, ScenTree & outTree,
                   CopulaDef::CopInfoForecastErrors const & tgCop,
                   DimT const nBr, DimT const dT = 1)
{
	if (initTree.is_empty()) {

	}
}
