/// definitions for scenario generation using forecast errors

#ifndef COP_GEN_FC_ERR_HPP
#define COP_GEN_FC_ERR_HPP

#include "copula-info.hpp"

namespace CopulaDef {

// ----------------------------------------------------------------------------
/// target copula described by historical data and forecasts
/**
	This just a version of \a CopInfoData with different constructors
	and (possibly) extra post-processing routines.
	The final result is a fan of values for a given number of periods ahead.

	Note that the variables we send to \a CopInfoData should be sorted
	by length of the forecast, so we can later extend this to multi-stage trees.
**/
class CopInfoForecastErrors : public CopInfoData {
private:
	DimT N;        ///< number of 'real' variables; note that nVar = N * T
	DimT T;        ///< number of stoch. periods (1..T; root/now = 0)
	MatrixD fCast; ///< forecast for the whole time horizon [nVar, T]

	/// override the default version from the base class
	/**
		This function is called from the constructor of the base class,
		which is too early. We therefore keep this empty and create
		a new version to be called from the constructors here.
	**/
	void setup_2d_targets() override {}

	/// creates objects for the 2D targets; called from the constructors
	/**
		\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (i, dt ± u) for u = 1..perVarDt
		\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (j, dt ± u) for u = 0..intVarDt
	**/
	void setup_2d_targets(int perVarDt, int intVarDt);

	inline DimT rowOf(DimT i, DimT dt) { return (dt - 1) * N + i; }

public:
	/// constructor with the historical forecast errors
	/**
		\param[in] forecast  forecast for the whole time horizon [N, T]
		\param[in] histFErr  historical forecast errors [N * T, nPts];
		                     columns are grouped by periods: val(i,dt) for
		                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
		\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (i, dt ± u) for u = 1..perVarDt
		\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (j, dt ± u) for u = 0..intVarDt
	**/
	CopInfoForecastErrors(MatrixD const & forecast,
	                      MatrixD const & histFErr,
	                      int perVarDt = 1, int intVarDt = 0);

	/// different types of sorting of the input historical data
	/**
		In all cases, the historical data includes values \c val(i,t) for
		variable \c i at time \c t and forecasts \c fc(i,t,t+dt) of value
		\c val(i,t+dt) made at time \c t.
		Row for time \c t always includes \c val(i,t), but may include
		different forecast. In addition, the data might be sorted ascending
		(old-to-new) or descending (new-to-old) in time.
	*/
	enum class HistDataSort {
		fCastTimeAsc,  ///< row(t): \c f(i,t,t+dt) for all dt; old-to-new
		fCastTimeDesc, ///< row(t): \c f(i,t,t+dt) for all dt; new-to-old
		valueTimeAsc,  ///< row(t): \c f(i,t-dt,t) for all dt; old-to-new
		valurTimeDesc  ///< row(t): \c f(i,t-dt,t) for all dt; new-to-old
	};

	/// constructor with the historical forecast errors
	/**
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
	**/
	CopInfoForecastErrors(MatrixD const & forecast,
	                      MatrixD const & histData,
	                      HistDataSort const dataSort,
	                      int perVarDt = 1, int intVarDt = 0);

	/// get the matrix of historical forecast errors
	/**
		This is needed if one uses the constructor with historical data
		and forecast, where the historical errors get computed.
	**/
	MatrixD const & hist_forecast_errors() const { return hData; }

	/// convert scenarios of errors to scenarios of the original values
	/**
		\param[in]  errSc  error-scenarios; [nVars, nSc]
		\param[out] scens  output scenarios; nSc * [T, N]
	**/
	void errors_to_values(MatrixD const & errSc,
	                      std::vector<MatrixD> & scens) const;
};

} // namespace CopulaDef


namespace FcErr_Gen {

} // namespace FcErr_Gen

#endif // COP_GEN_FC_ERR_HPP
