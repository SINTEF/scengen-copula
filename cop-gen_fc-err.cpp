#include "cop-gen_fc-err.hpp"

//#include <iostream>
#include <fstream>
#include <cstdlib>  // for system() command
#include <cstdio>   // for remove() command
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


// return the number of forecast steps in historical data
/*
	\param[in]  histData  historical values and forecasts [nPts, nVars];
	\param[in]   dataFmt  data format of \c histData
	\param[in]         N  the original dimension (number of orig. variables)
*/
DimT FcErr_Gen::nmb_dts_in_data(MatrixD const & histData,
                                HistDataFormat const dataFmt, DimT const N)
{
	DimT T;
	if ((dataFmt & HistDataFormat::hasErrors) == HistDataFormat::hasErrors) {
		// data includes errors -> one column per variable and time step
		T = histData.size2() / N;
		if (histData.size2() != N * T)
			throw std::length_error
				("Inconsistent dimensions of input matrices.");
	} else {
		// forecasts and values -> 1 extra column per variable
		T = histData.size2() / N - 1;
		if (histData.size2() != N * (T + 1))
			throw std::length_error
				("Inconsistent dimensions of input matrices.");
	}
	return T;
}


/*
	\param[in]  histData  historical values and forecasts [nPts, nPts];
	\param[out] histFErr  the computed historical forecast errors [nVars, nPts]
	\param[in]         N  the original dimension (number of orig. variables)
	\param[in]   dataFmt  data format of \c histData
	\param[in]  errTypes  types of errors, per variable [nVars]
	                      empty list means linear differences
*/
void FcErr_Gen::process_hist_data(MatrixD const & histData, MatrixD & histFErr,
                                  DimT const N, HistDataFormat const dataFmt,
                                  FcErrTypeList const & errTypes)
{
	HistDataFormat curFmt = dataFmt; // current format, updated on the way

	if ((dataFmt & HistDataFormat::newToOld) == HistDataFormat::newToOld) {
		throw std::logic_error
			("new-to-old sorting of historical data is not supported yet");
	}
	if ((dataFmt & HistDataFormat::rowPerVal) == HistDataFormat::rowPerVal) {
		throw std::logic_error("row-per-value format is not supported yet");
	}

	DimT T;     // number of periods
	DimT nPts;  // number of observations in the data
	DimT nVars; // number of generated variables; = N * T
	DimT i, dt, cF, cV, rE, j;

	if ((curFmt & HistDataFormat::hasErrors) == HistDataFormat::hasErrors) {
		// histData includes errors
		T = histData.size2() / N;
		nPts = histData.size1();
		nVars = N * T;
		histFErr.resize(nVars, nPts, false); // transposed rel. to histData!

		if ((curFmt & HistDataFormat::grByPer) == HistDataFormat::grByPer) {
			// correct grouping -> just copy the data
			for (i = 0; i < nVars; ++i)
				for (j = 0; j < nPts; ++j)
					histFErr(i, j) = histData(j, i);
		} else {
			// need to reorder the columns
			for (dt = 0; dt < T; ++dt) {
				for (i = 0; i < N; ++i) {
					rE = dt * N + i; // var. row in histFErr
					cV = i * T + dt; // var. column in histData
					for (j = 0; j < nPts; ++j) {
						histFErr(rE, j) = histData(j, cV);
					}
				}
			}
			curFmt |= HistDataFormat::grByPer; // update the indicator
		}
	}
	else
	{
		// the data are given as values + forecasts
		T = histData.size2() / N - 1;
		nPts = histData.size1() - T; // need extra data to compute the errors
		nVars = N * T;
		histFErr.resize(nVars, nPts, false); // transport rel. to histData!

		// local version of errTypes, allowing default value if empty
		FcErrTypeList const & fcErrTypes
			= errTypes.size() > 0 ? errTypes
			                      : FcErrTypeList(N, FcErrorType::linDiff);

		for (i = 0; i < N; ++i) {
			cV = i * (T+1); // column of the variable value
			for (dt = 1; dt <= T; ++dt) {
				// NB: this assumes the default input-file format!
				cF = cV + dt;      // column of the forecast
				rE = hist_data_row_of(i, dt, N, T); // row within histFErr
				switch (fcErrTypes[i]) {
				case FcErrorType::linDiff:
					for (j = 0; j < nPts; ++j) {
						// error = value - forecast
						histFErr(rE, j) = histData(j + dt, cV) - histData(j, cF);
					}
					break;
				case FcErrorType::relDiff:
					for (j = 0; j < nPts; ++j) {
						// error = (value-forecast) / forecast = value/forecast - 1
						histFErr(rE, j) = histData(j + dt, cV) / histData(j, cF) - 1;
					}
					break;
				case FcErrorType::revRelDiff:
					for (j = 0; j < nPts; ++j) {
						// error = (forecast - value) / value = forecast/value - 1
						histFErr(rE, j) = histData(j, cF) / histData(j + dt, cV) - 1;
					}
					break;
				default:
					throw std::logic_error("unknown type of forecast error");
				}
			}
		}
		curFmt |= HistDataFormat::hasErrors; // update the indicator
		curFmt |= HistDataFormat::grByPer;   // update the indicator
	}

	assert (curFmt == (HistDataFormat::grByPer | HistDataFormat::hasErrors)
	        && "check that we have the correct format at the end");
}


// convert scenarios of errors to scenarios of the original values
/*
	\param[in]  errSc    error-scenarios; [nVars, nSc]
	\param[in] forecast  forecast for the whole time horizon; [T, N]
	\param[out] scens    output scenarios; nSc * [T, N]
	\param[in] errTypes  types of errors, per variable [nVars]
	                     empty list means linear differences
*/
void FcErr_Gen::errors_to_scens(MatrixD const & errSc, MatrixD const & forecast,
                                 std::vector<MatrixD> & scens,
                                 FcErrTypeList const & errTypes)
{
	assert (errSc.size1() * errSc.size2() > 0 && "check non-empty input");
	assert (forecast.size1() * forecast.size2() > 0 && "check non-empty input");
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

	// local version of errTypes, allowing default value if empty
	FcErrTypeList const & fcErrTypes
		= errTypes.size() > 0 ? errTypes
		                      : FcErrTypeList(N, FcErrorType::linDiff);

	scens.resize(nSc);
	for (DimT s = 0; s < nSc; ++s) {
		scens[s].resize(T, N);
		for (DimT t = 0; t < T; ++t) {
			for (DimT i = 0; i < N; ++i) {
				DimT r = N * t + i; // row in errSc
				switch (fcErrTypes[i]) {
				case FcErrorType::linDiff:
					scens[s](t, i) = forecast(t, i) + errSc(r, s);
					break;
				case FcErrorType::relDiff:
					scens[s](t, i) = forecast(t, i) * (1 + errSc(r, s));
					break;
				case FcErrorType::revRelDiff:
					scens[s](t, i) = forecast(t, i) / (1 + errSc(r, s));
					break;
				default:
					throw std::logic_error("unknown type of forecast error");
				}
			}
		}
	}
}

// convert generated errors (2-stage) to scenarios of the error terms
/*
	\param[in]    errSc  error-scenarios; [nVars, nSc]
	\param[in]        N  number of the original variables
	\param[out]   scens  output scenarios; nSc * [T, N]
*/
void FcErr_Gen::errors_to_scens(MatrixD const & errSc, DimT const N,
                                std::vector<MatrixD> & scens)
{
	assert (errSc.size1() * errSc.size2() > 0 && "check non-empty input");
	DimT nVars = errSc.size1();
	DimT nSc = errSc.size2();
	DimT T = nVars / N;
	if (nVars != T * N)
		throw std::length_error
			("inconsistent sizes of error and forecast matrices");

	DimT i, r, s, t;

	scens.resize(nSc);
	for (s = 0; s < nSc; ++s) {
		scens[s].resize(T, N);
		for (t = 0; t < T; ++t) {
			for (i = 0; i < N; ++i) {
				r = N * t + i; // row in errSc
				scens[s](t, i) = errSc(r, s);
			}
		}
	}
}


// --------------------------------------------------------------------------
// class CopInfoForecastErrors

// the 'common constructor', used for delegation (never called by itself)
/*
	\param[in]  nmbVars  number of variables
	\param[in] histFErr  historical forecast errors [N * T, nPts];
	                     columns are grouped by periods: val(i,dt) for
	                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
*/
CopInfoForecastErrors::CopInfoForecastErrors(DimT const nmbVars,
	                                         MatrixD const & histFErr)
: CopInfoData(histFErr, false), // 'false' prohibits call of setup_2d_targets()
  N(nmbVars), T(histFErr.size1() / N)
{
	if (histFErr.size1() != N * T)
		throw std::length_error("inconsistent dimensions of input data");
}

// constructor with the list of variable pairs for the copulas
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
CopInfoForecastErrors::CopInfoForecastErrors(
	DimT const nmbVars, MatrixD const & histFErr,
	std::list<std::pair<DimT, DimT>> const & copVars)
: CopInfoForecastErrors(nmbVars, histFErr)
{
	CopInfoData::setup_2d_targets(copVars);
}

// setup target 2d copulas, based on a matrix of max dts.
void CopInfoForecastErrors::setup_2d_targets(ublas::symmetric_matrix<int> const & maxDts)
{
	// create the list of copulas:
	std::list<std::pair<DimT, DimT>> copVars;

	// NB: the first two loops were reversed in prev. version -> CHECK!
	for (DimT dt1 = 1; dt1 <= T; ++ dt1) {
		for (DimT i = 0; i < N; ++i) {
			DimT v1 = row_of(i, dt1);
			for (DimT j = 0; j < N; ++j) {
				for (int dt2 = dt1 - maxDts(i, j); dt2 <= (int) dt1 + maxDts(i, j); ++dt2) {
					if (dt2 >= 1 && dt2 <= (int) T) {
						DimT v2 = row_of(j, dt2);
						if (v2 > v1) {
							MSG(TrInfo2, "adding 2d-copula for variables ("
							    << i << "," << dt1 << ") and (" << j << "," << dt2 << ")");
							copVars.push_back(std::make_pair(v1, v2));
						}
					}
				}
			}
		}
	}

	CopInfoData::setup_2d_targets(copVars);
}

// constructor with matrix of max-dt values
/*
	\param[in]  nmbVars  number of variables
	\param[in] histFErr  historical forecast errors [N * T, nPts];
	                     columns are grouped by periods: val(i,dt) for
	                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
	\param[in]   maxDts  matrix of max-dt for all variable pairs
*/
CopInfoForecastErrors::CopInfoForecastErrors(DimT const nmbVars,
	                                         MatrixD const & histFErr,
	                                         ublas::symmetric_matrix<int> const & maxDts)
: CopInfoForecastErrors(nmbVars, histFErr)
{
	setup_2d_targets(maxDts);
}

// constructor with max perVarDt and intVarDt for 2D-copula creation
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
: CopInfoForecastErrors(nmbVars, histFErr)
{
	TRACE(TrDetail, "CopInfoForecastErrors constructor with perVarDt and intVarDt");

	// setup the matrix of dts
	ublas::symmetric_matrix<int> maxDts(N);

	for (DimT i = 0; i < N; ++i) {
		maxDts(i, i) = perVarDt;
		for (DimT j = 0; j < i; ++j)
			maxDts(i, j) = intVarDt;
	}
	DBGSHOW_NL(TrInfo2, maxDts);

	setup_2d_targets(maxDts);
}


// version of \c hist_data_row_of() for use inside the class (fewer parameters)
DimT CopInfoForecastErrors::row_of(DimT const i, DimT const dt)
{
	return hist_data_row_of(i, dt, N, T);
}


// --------------------------------------------------------------------------
// class ScenTree

void ScenTree::fill_vals_per_var()
{
	if (scenVals.size() == 0) {
		std::cerr << "called fill_vals_per_var without output data - ignoring"
		          << std::endl;
		return;
	}
	DimT s, S = scenVals.size();
	DimT t, T = scenVals[0].size1();
	DimT i, N = scenVals[0].size2();
	DimT r, R = T + (rootVals.size() == N ? 1 : 0);
	varVals.resize(N);
	for (i = 0; i < N; ++i) {
		varVals[i].resize(R, S);
		r = 0;
		if (rootVals.size() == N) {
			for (s = 0; s < S; ++s)
				varVals[i](r, s) = rootVals(i);
			++r;
		}
		for (t = 0; t < T; ++t)
			for (s = 0; s < S; ++s)
				varVals[i](t + r, s) = scenVals[s](t, i);
	}
}


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

void ScenTree::set_root_values(std::vector<double> const & curVal)
{
	DimT N = curVal.size();
	if (scenVals.size() > 0 && N != scenVals[0].size2())
		throw std::length_error("vector of root values has a wrong length");
	DimT i;
	rootVals.resize(N);
	for (i = 0; i < N; ++i)
		rootVals(i) = curVal[i];
}


// add a fixed scenario to a generated tree
std::pair<DimT, DimT> ScenTree::add_scenario(MatrixD const & scVals, double const prob)
{
	DimT S = scenVals.size();
	if (S > 0) {
		// values from before - check dimensions
		if (scVals.size1() != scenVals[0].size1()
			|| (scVals.size2() != scenVals[0].size2()))
		throw std::logic_error("Wrong dimension of the new scenarios.");
	}
	scenVals.push_back(scVals);

	// (re)calculate probabilities
	if (scProbs.size() == S) {
		for (auto it = scProbs.begin(); it != scProbs.end(); ++it)
			*it *= (1 - prob);
		scProbs.resize(S + 1, true);
	} else {
		if (scProbs.size() != 0)
			throw std::logic_error("Non-zero vector scProb has wrong size.");
		scProbs = VectorD(S + 1, (1.0 - prob) / S);
	}
	scProbs[S] = prob;

	if (varVals.size() > 0) {
		// re-generate the per-value description
		fill_vals_per_var();
	}

	// find which scenario the new one branches off and when
	std::pair<DimT, DimT> branchOff = std::make_pair(0, 0);
	DimT T = nmb_periods();
	DimT N = nmb_variables();
	for (DimT sc = 0; sc < S; ++sc) {
		DimT i, t;  // they are both used outside their for-loops
		for (t = 0; t < T; ++t) {
			for (i = 0; i < N; ++i) {
				if (fabs(scenVals[sc](t, i) - scVals(t, i)) > 1e-6)
					break;
			}
			if (i < N) {
				// difference in values
				break;
			}
		}
		// now, t is the first time where the new scenario differs from sc
		if (t > branchOff.second) {
			branchOff.first = sc;
			if (branchOff.second < T) {
				branchOff.second = t;
			} else {
				WARNING("new scenario is equal to an existing scen. " << sc);
				branchOff.second = T-1;  // alt: merge the two scenarios
			}
		}
		if (t >= T-1) {
			// it cannot be bigger -> stop the search
			break;
		}
	}

	return branchOff;
}


// get values for a given scenarios; matrix [T, N]
MatrixD const & ScenTree::values_for_scen(DimT const sc) const
{
	return scenVals[sc];
}

// get values for a given variable; matrix [T, S]
MatrixD const & ScenTree::values_for_var(DimT const i)
{
	if (varVals.size() == 0)
		fill_vals_per_var();
	return varVals[i];
}

/// get scenario probabilities
VectorD const & ScenTree::scen_probs() const
{
	if (scProbs.size() > 0 && scProbs.size() != scenVals.size()) {
		throw std::logic_error("Non-zero vector scProb has wrong size.");
	}
	return scProbs;
}


void ScenTree::display_per_scen(std::ostream & outStr) const
{
	DimT s, S = scenVals.size();
	for (s = 0; s < S; ++s) {
		outStr << "# values for scenario " << s+1 << ":" << std::endl;
		if (rootVals.size() == scenVals[0].size2())
			outStr << rootVals << std::endl;
		outStr << scenVals[s] << std::endl;
	}
}

void ScenTree::display_per_var(std::ostream & outStr,
                               bool const useVarHeader,
                               std::string const varSep)
{
	if (varVals.size() == 0)
		fill_vals_per_var();

	DimT i, N = scenVals[0].size2();
	for (i = 0; i < N; ++i) {
		if (i > 0)
			outStr << varSep;
		if (useVarHeader)
			outStr << "# scenarios for variable " << i+1 << ":" << std::endl;
		outStr << varVals[i];
	}
}

int ScenTree::make_gnuplot_charts(std::string const & baseFName,
                                  MatrixD const * const p2forecast,
                                  std::string const & gnuplotExe)
{
	if (gnuplotExe.size() == 0) {
		WARNING("called make_gnuplot_charts() with empty gnuplotExe");
		return 1;
	}

	DimT s, S = scenVals.size();
	DimT    T = scenVals[0].size1();
	DimT i, N = scenVals[0].size2();
	bool hasFcast = p2forecast && (p2forecast->size1()*p2forecast->size2() > 0);

	// file for the script
	std::string scrFName = "make_plots.gp";
	/*{
		std::stringstream scrFNameStr;
		scrFNameStr << "make_plots_" << time(nullptr) << ".gp";
		scrFName = scrFNameStr.str();
	}*/
	std::ofstream scrF(scrFName, std::ios::out);
	if (!scrF)
		throw std::ios_base::failure("could not open output file " + scrFName);

	// file for the output scenario data
	std::string outFName = "out_scens.out";
	/*{
		std::stringstream outFNameStr;
		outFNameStr << "out_scens_" << time(nullptr) << ".out";
		outFName = outFNameStr.str();
	}*/
	std::ofstream outF(outFName, std::ios::out);
	if (!outF)
		throw std::ios_base::failure("could not open output file " + outFName);

	// first, write the data file
	if (hasFcast) {
		outF << "# forecast" << '\n';
		if (rootVals.size() == N)
			outF << rootVals << '\n';
		outF << *p2forecast << "\n\n";
	}
	display_per_var(outF, true, "\n\n");
	outF.close();

	// now write the script;
	scrF << "set xrange [-0.1:" << T + 0.1 << "]" << '\n';
	scrF << "set xtics 1" << '\n';
	DimT figW = 800;
	if (S > 24)
		scrF << "unset key" << '\n';  // too many scenarios to show a legend
	else {
		scrF << "set key outside right" << '\n';
		figW += 160;
	}
	scrF << "set term png enhanced size " << figW << ",600 truecolor" << '\n';
	scrF << '\n';
	for (i = 0; i < N; ++i) {
		string figFName
			= baseFName + (varNames.size() == N ? varNames[i]
			                                    : std::to_string(i+1)) + ".png";
		scrF << "set output '" << figFName << "'" << '\n';
		scrF << "set title 'Scenarios for variable " << i+1 << "'" << '\n';

		scrF << "plot \\" << '\n';
		//
		std::string xCol = ((rootVals.size() == N ? "0" : "($0+1)"));
		if (hasFcast) {
			for (s = 0; s < S; ++s)
				scrF << "\t'" << outFName << "' using " << xCol << ":" << s+1
				     << " index " << i+1 << " w linespoints title 'scen " << s+1
				     << "', \\" << '\n';
			scrF << "\t'" << outFName << "' using " << xCol << ":" << i+1
			     << " index 0 w linespoints ls -1 lw 2 title 'forecast'\n";
		} else {
			for (s = 0; s < S; ++s)
				scrF << "\t'" << outFName << "' using " << xCol << ":" << s+1
				     << " index " << i << " w linespoints title 'scen " << s+1
				     << "'" << (s < S-1 ? ", \\" : "") << '\n';
		}
		//scrF << "pause -1 'press any key'\n";
		scrF << "set output" << '\n';
		scrF << '\n';
	}
	scrF.close();

	std::stringstream gpCmd;
	gpCmd << gnuplotExe << " " << scrFName;
	if (system(nullptr) == 0) {
		std::cerr << "Error: No command-line processor available "
		             "-> cannot execute gnuplot.\n";
	} else {
		TRACE(TrDetail, "Executing system command '" << gpCmd.str() << "'");
		int retV = system(gpCmd.str().c_str());
		if (retV != 0) {
			WARNING("Gnuplot command not successful, error code " << retV);
		}
	}

	// delete the temp. files
	remove(outFName.c_str());
	remove(scrFName.c_str());

	return 0;
}


// --------------------------------------------------------------------------
// class FcErrTreeGen

// minimal constructor, used only for delegating
/*
	\param[in]     nmbVars  number of stoch. variables
	\param[in]    histData  historical values and forecasts [nPts, nVars];
	\param[in]  dataFormat  sorting/structure of \c histData
	\param[in]    errTypes  list of error types [nVars]
	\param[in]   zeroMeans  enforce zero means for all errors?
*/
FcErrTreeGen::FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
                           HistDataFormat const dataFormat,
                           FcErrTypeList const & errTypes,
                           bool const zeroMeans)
: N(nmbVars), fErrTypes(errTypes), zeroErrorMeans(zeroMeans)
{
	if ((dataFormat & HistDataFormat::newToOld) == HistDataFormat::newToOld) {
		throw std::logic_error
			("new-to-old sorting of historical data is not supported yet");
	}
	if ((dataFormat & HistDataFormat::hasErrors) == HistDataFormat::hasErrors) {
		// histData includes historical forecast errors
		T = histData.size2() / N;
		if (histData.size2() != N * T)
			throw std::length_error("inconsistent dimensions of input data");
		_histFErr = histData; // copy the data
	} else {
		// histData includes observed value and forecasts
		T = histData.size2() / N - 1;
		if (histData.size2() != N * (T+1))
			throw std::length_error("inconsistent dimensions of input data");
	}

	// convert to the desired format: HistDataFormat(::grByPer & ::hasErrors)
	// note: We cannot use histFErr here, since it is a ref. to const matrix.
	//       However, it works with _histFErr, which is a non-const matrix
	process_hist_data(histData, _histFErr, N, dataFormat, errTypes);

	DBGSHOW_NL(TrDetail, histFErr);
}

// constructor one error type (same for all variables) and max*VarDt
/*
	\param[in]     nmbVars  number of stoch. variables
	\param[in]    histData  historical values and forecasts [nPts, nVars];
	\param[in]    dataSort  sorting/structure of \c histData
	\param[in]     errType  error type, used for all variables
	\param[in] maxPerVarDt  for var. (i, dt), we generate 2D-copulas with
	                        (i, dt ± u) for u = 1..perVarDt
	\param[in] maxIntVarDt  for var. (i, dt), we generate 2D-copulas with
	                        (j, dt ± u) for u = 0..intVarDt
	\param[in]   zeroMeans  enforce zero means for all errors?
*/
FcErrTreeGen::FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
                           HistDataFormat const dataFormat,
                           FcErrorType const errType,
                           int const maxPerVarDt, int const maxIntVarDt,
                           bool const zeroMeans)
: FcErrTreeGen(nmbVars, histData, dataFormat, FcErrTypeList(nmbVars, errType), zeroMeans)
{
	perVarDt = maxPerVarDt;
	intVarDt = maxIntVarDt;
}

// constructor variable error types and specified copula pairs
/*
	\param[in]     nmbVars  number of stoch. variables
	\param[in]    histData  historical values and forecasts [nPts, nVars];
	\param[in]    dataSort  sorting/structure of \c histData
	\param[in]    errTypes  list of error types [nVars]
	\param[in]      maxDts  matrix of max-dt for all variable pairs
	\param[in]   zeroMeans  enforce zero means for all errors?
*/
FcErrTreeGen::FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
                           HistDataFormat const dataFormat,
                           FcErrTypeList const & errorTypes,
                           ublas::symmetric_matrix<int> const & maxDts,
                           bool const zeroMeans)
: FcErrTreeGen(nmbVars, histData, dataFormat, errorTypes, zeroMeans)
{
	copMaxDt = maxDts;
}


// generate a 2-stage tree (a fan), with given number of scenarios
/*
	\param[in]    nScen  number of scenarios to generate
	\param[out] outTree  the resulting scenario tree
*/
void FcErrTreeGen::gen_2stage_tree(DimT const nScen, ScenTree & outTree)
{
	assert (N > 0 && "sanity check");

	outLvl = (OutputLevel) ((int) outLvl + 1);
	CopInfoBy2D::Ptr p2tgCop;
	if (copVarPairs.size() > 0) {
		p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, histFErr, copVarPairs);
	} else {
		if (copMaxDt.size1() > 0)
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, histFErr, copMaxDt);
		else
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, histFErr, perVarDt, intVarDt);
	}
	outLvl = (OutputLevel) ((int) outLvl - 1);

	auto p2tgMargs = boost::make_shared<MarginDistrib::SampleMargins>
		(histFErr);
	if (zeroErrorMeans)
		p2tgMargs->fixMeans(0.0);

	CopulaScen::CopulaSample copSc(p2tgCop, nScen);
	copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);

	// transform margins to the target distribution
	MatrixD errScens;  // dim = [nVars, nSc]
	p2tgMargs->get_margin_distr(copRanks, errScens);
	//DISPLAY_NL(errScens);

	// post-process the error scenarios
	if (errScenPP)
		errScenPP->process(errScens);

	// convert to the scenarios for the original multi-period problem
	if (curFcast.get().size1() > 0)
		errors_to_scens(errScens, curFcast, outTree.scenVals, fErrTypes);
	else
		errors_to_scens(errScens, N, outTree.scenVals);
}


// filter copVarPairs to include only variables up to tMax
std::list<std::pair<DimT, DimT>> FcErrTreeGen::cop_vars_sublist(DimT tMax)
{
	std::list<std::pair<DimT, DimT>> sublist;
	DimT iMax = N * tMax;

	for (auto&& cop : copVarPairs) {
		if (cop.first < iMax && cop.second < iMax)
			sublist.push_back(cop);
	}

	return sublist;
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
	MSG(TrDetail2, "historical data for the new stage:\n" << errM);

	outLvl = (OutputLevel) ((int) outLvl + 1);
	CopInfoForecastErrors::Ptr p2tgCop;
	if (copVarPairs.size() > 0) {
		if (nPers < T) {
			// have to filter the list of variables
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(
				N, errM, cop_vars_sublist(nPers));
		} else {
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, errM, copVarPairs);
		}
	} else {
		if (copMaxDt.size1() > 0)
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, errM, copMaxDt);
		else
			p2tgCop = boost::make_shared<CopInfoForecastErrors>(N, errM, perVarDt, intVarDt);
	}
	outLvl = (OutputLevel) ((int) outLvl - 1);

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
	if (zeroErrorMeans)
		p2tgMargs->fixMeans(0.0);
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
		errors_to_scens(errScens, *p2forecast, p2outTree->scenVals);
	}
	*/
}


// generate a regular multistage tree, with given branching per period
/*
	\param[in] branching  branching factor per period
	\param[out]  outTree  the resulting scenario tree

	\note By a regular tree we mean a scenario tree where all nodes
	      in the same period have the same number of child nodes.
*/
void FcErrTreeGen::gen_reg_tree(VectorI const & branching, ScenTree & outTree)
{
	DimT nmbPer = T;  // number of periods we generate scenarios for
	if (curFcast.get().size1() > 0) {  // NB: curFcast is a reference_wrapper
		if (curFcast.get().size2() != N)
			throw std::length_error("wrong number of variables in the forecast");
		if (branching.size() != curFcast.get().size1())
			WARNING("different number of periods in branching vector and forecast");
		if (curFcast.get().size1() < nmbPer)
			nmbPer = curFcast.get().size1();
	}
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

	// post-process the error scenarios
	if (errScenPP)
		errScenPP->process(errScens);

	if (curFcast.get().size1() > 0) {
		DBGSHOW_NL(TrDetail, curFcast.get());
		errors_to_scens(errScens, curFcast, outTree.scenVals, fErrTypes);
	} else {
		errors_to_scens(errScens, N, outTree.scenVals);
	}
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
