#include "copula-info.hpp"

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
: CopInfoData(histFErr), T(forecast.size2()), fCast(forecast)
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
	nPts = histData.size1() - T; // need extra data to compute the errors
	DimT i, dt, cF, cV, rE, j;

	// fill the hData matrix (NB: transposed rel. to histData!)
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

	// rest of the setup
	fill_ranks_etc();                     // fill hRanks and hU01
	setup_2d_targets(perVarDt, intVarDt); // create the bivariate copula objects
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
