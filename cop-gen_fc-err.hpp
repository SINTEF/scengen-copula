/// definitions for scenario generation using forecast errors

#ifndef COP_GEN_FC_ERR_HPP
#define COP_GEN_FC_ERR_HPP

#include "copula-info.hpp"
#include "margins.hpp"
#include "copula-sample.hpp"

#include "bitmask_operators.hpp"

#include <iostream>


/// classes and methods specific for the forecast-error-based generator
namespace FcErr_Gen {

/// different formats of historical input data
/**
	This includes binary indicators for format of the historical data.
	These can be combined used the templates from bitmask_operators.hpp.

	In all cases, the data have variables in columns and observations in rows.

	If ::hasErrors==0, then the data consists of the actual observed values,
	followed by forecasts. Timing of the forecasts is controlled by ::rowPerVal.
**/
enum class HistDataFormat {
	standard  = 0, ///< standard/default format
	grByPer   = 1, ///< values grouped by period; default = grouped by var.
	hasErrors = 2, ///< data consists of errors; default = values and forecasts
	newToOld  = 4, ///< values are sorted new-to-old; default = old-to-new
	rowPerVal = 8  ///< row(t) includes forecasts for value at time t;
	               ///< default = row(t) includes forecasts made at time t;
	               ///< this is valid only if hasErrors = 0
};

/// different ways of computing forecast errors
/**
	The two types of relative differences can be useful if we know
	that either forecasts or values are never zero.
**/
enum class FcErrorType {
	linDiff,    ///< linear/normal difference: E = X - F
	relDiff,    ///< relative difference: E = (X - F) / F
	revRelDiff  ///< reversed rel. difference: E = (F - X) / X
};
using FcErrTypeList = std::vector<FcErrorType>;


/// base class for post-processing of scenarios
class ScenPostProcess {
public:
	virtual void process(MatrixD & scens) const = 0;  // dim = [nVars, nSc])
};

// ----------------------------------------------------------------------------
/// target copula described by historical data and forecasts
/**
	This just a version of \a CopInfoData with different constructors
	and (possibly) extra post-processing routines.
	The final result is a fan of values for a given number of periods ahead.

	Note that the variables we send to \a CopInfoData should be sorted
	by length of the forecast, so we can later extend this to multi-stage trees.

	\todo If the forecast length varies a lot, we might want to use more recent
	      data for the shortest forecasts. To achieve this, we would need to
	      have another class derived from CopInfoBy2D with bivariate 2D
	      copulas - this would allow us to have different lengths.
	      For the input file, the matrix-reading code now handles 'NaN'
	      in the input, which is converted to signaling_NaN. All we need to do
	      is to check for NaN when reading the matrix.
**/
class CopInfoForecastErrors : public CopulaDef::CopInfoData {
private:
	DimT N;        ///< number of 'real' variables; note that nVars = N * T
	DimT T;        ///< number of stoch. periods (1..T; root/now = 0)
	//MatrixD fCast; ///< forecast for the whole time horizon [nVar, T]

	/// version of \c hist_data_row_of() for use inside the class;
	/// fewer parameters, just to save typing
	inline DimT row_of(DimT const i, DimT const dt);

	/// constructor with the common parts; used only for delegating
	/**
		\param[in]  nmbVars  number of variables
		\param[in] histFErr  historical forecast errors [N * T, nPts];
		                     columns are grouped by periods: val(i,dt) for
		                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
	**/
	CopInfoForecastErrors(DimT const nmbVars, MatrixD const & histFErr);

	//using CopulaDef::CopInfoData::setup_2d_targets;

	/// setup target 2d copulas, based on a matrix of max dts.
	void setup_2d_targets(ublas::symmetric_matrix<int> const & maxDts);

public:
	/// constructor with perVarDt and intVarDt
	/**
		\param[in]  nmbVars  number of variables
		\param[in] histFErr  historical forecast errors [N * T, nPts];
		                     columns are grouped by periods: val(i,dt) for
		                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
		\param[in] perVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (i, dt ± u) for u = 1..perVarDt
		\param[in] intVarDt  for var. (i, dt), we generate 2D-copulas with
		                     (j, dt ± u) for u = 0..intVarDt
	**/
	CopInfoForecastErrors(DimT const nmbVars, MatrixD const & histFErr,
	                      int perVarDt, int intVarDt);

	/// constructor with matrix of max-dt values
	/**
		\param[in]  nmbVars  number of variables
		\param[in] histFErr  historical forecast errors [N * T, nPts];
		                     columns are grouped by periods: val(i,dt) for
		                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
		\param[in]   maxDts  matrix of max-dt for all variable pairs
	**/
	CopInfoForecastErrors(DimT const nmbVars, MatrixD const & histFErr,
	                      ublas::symmetric_matrix<int> const & maxDts);

	/// constructor with a list of pairs to create
	/**
		\param[in]  nmbVars  number of variables
		\param[in] histFErr  historical forecast errors [N * T, nPts];
		                     columns are grouped by periods: val(i,dt) for
		                     (0,1), (1,1),...,(N,1), (0,2),...,(N-1,T)
		\param[in]  copVars  list of pairs of variables to create a copula for
	**/
	CopInfoForecastErrors(DimT const nmbVars, MatrixD const & histFErr,
	                      std::list<std::pair<DimT, DimT>> const & copVars);


	/// get the matrix of historical forecast errors
	/**
		This is needed if one uses the constructor with historical data
		and forecast, where the historical errors get computed.
	**/
	MatrixD const & hist_forecast_errors() const { return hData; }
};


/// return the number of forecast steps in historical data
/**
	\param[in]  histData  historical values and forecasts [nPts, nVars];
	\param[in]   dataFmt  data format of \c histData
	\param[in]         N  the original dimension (number of orig. variables)
**/
DimT nmb_dts_in_data(MatrixD const & histData, HistDataFormat const dataFmt,
                     DimT const N);

/// process historical data (compute historical forecast errors etc)
/**
	\param[in]  histData  historical values and forecasts [nPts, nVars];
	\param[out] histFErr  the computed historical forecast errors [nVars, nPts]
	\param[in]         N  the original dimension (number of orig. variables)
	\param[in]   dataFmt  data format of \c histData
	\param[in]  errTypes  types of errors, per variable [nVars];
	                      empty list means linear differences

	The output of \c histFErr is transposed ([nVars, nPts]) and its format
	should be \c (HistDataFormat::grByPer & HistDataFormat::::hasErrors).
**/
void process_hist_data(MatrixD const & histData, MatrixD & histFErr,
                       DimT const N, HistDataFormat const dataFmt,
                       FcErrTypeList const & errTypes = FcErrTypeList());


/// convert scenarios of errors to scenarios of the original values
/**
	\param[in]    errSc  error-scenarios; [nVars, nSc]
	\param[in] forecast  forecast for the whole time horizon; [T, N]
	\param[out]   scens  output scenarios; nSc * [T, N]
	\param[in] errTypes  types of errors, per variable [nVars]
	                     empty list means linear differences

	This version takes the generated 1-period sample of prediction errors
	and combines them with the supplied forecast to create multi-period
	scenarios for the original values.
**/
void errors_to_scens(MatrixD const & errSc, MatrixD const & forecast,
                     std::vector<MatrixD> & scens,
                     FcErrTypeList const & errTypes = FcErrTypeList());

/// convert generated errors (2-stage) to scenarios of the error terms
/**
	\param[in]    errSc  error-scenarios; [nVars, nSc]
	\param[in]        N  number of the original variables
	\param[out]   scens  output scenarios; nSc * [T, N]

	This version simply re-arranges the values from the generated
	1-period samples to multi-period scenarios.
**/
void errors_to_scens(MatrixD const & errSc, DimT const N,
                     std::vector<MatrixD> & scens);


class FcErrTreeGen; // forward declaration


/// scenario tree for use in the forecast-error-based generator
class ScenTree {
friend class FcErrTreeGen;
private:
	/*
	DimT N; ///< dimension, i.e., the number of variables/margins
	DimT S; ///< number of scenarios
	DimT T; ///< number of periods
	*/

	/// \name structures for storing the tree values
	/**
		We store the scenario tree as a fan, with matrix of values per scenario.
		In addition, we add the possibility to have a value for common root
		and a branching information per period.
	**/
	///@{
		std::vector<MatrixD> scenVals; ///< scenario values; S * [T, N]
		VectorI branching; ///< optional branching per period; [T]
		VectorD rootVals;  ///< optional root values; [N]
		VectorD scProbs;   ///< optional scenario probabilities [S]

		std::vector<MatrixD> varVals; ///< values per var, for output; N * [T, S]
		void fill_vals_per_var();     ///< fill \c varVals from \c scenVals

		std::vector<string> varNames; ///< optional variable names (for output)
	///@}
public:
	ScenTree() = default;

	/// check if the tree is empty
	bool is_empty() const;

	/// get the number of scenarios
	DimT nmb_scens() const { return scenVals.size(); }

	/// get the number of periods
	DimT nmb_periods() const;

	/// get the number of variables
	DimT nmb_variables() const;

	/// set the current values
	void set_root_values(VectorD const & curVal) { rootVals = curVal; }
	/// set the current values
	void set_root_values(std::vector<double> const & curVal);

	/// add a fixed scenario to a generated tree
	/**
		- used to add a forecast to the error-based tree
		- must be called after scenario-generation
		- NB: the result is no longer a balanced tree!

		\param[in] scVals  matrix of scenario values; matrix [T, N]
		\return    pair <sc, t> denoting that the new scenario branches
		           off scen. sc in per t (i.e., is equal to sc up to t-1)
	**/
	std::pair<DimT, DimT> add_scenario(MatrixD const & scVals, double const prob);

	/// get values for a given scenarios; matrix [T, N]
	MatrixD const & values_for_scen(DimT const sc) const;

	/// get values for a given variable; matrix [T, S]
	/**
		NB: cannot be const, since it fills varVals if required!
	**/
	MatrixD const & values_for_var(DimT const i);

	/// get scenario probabilities; can return empty vector if equiprobable
	VectorD const & scen_probs() const;

	/// set variable names (for output)
	void set_var_names(std::vector<string> const & varNames)
	{ this->varNames = varNames; }

	/// simple stream printout of the values, one matrix per scenario
	void display_per_scen(std::ostream & outStr = std::cout) const;

	/// simple stream printout of the values, one matrix per variable
	/**
		\note Cannot be \c const, because it fills \c varVals if necessary
	**/
	void display_per_var(std::ostream & outStr = std::cout,
	                     bool const useVarHeader = true,
	                     std::string const varSep = "\n");

	/// make gnuplot charts, per variable
	/**
		\note Cannot be \c const, because it fills \c varVals if necessary

		At the moment, this simply calls gnuplot using the \c system() command.
		An alternative could be using the gnuplot-iostream interface from
		http://www.stahlke.org/dan/gnuplot-iostream/, but that would introduce
		extra dependency...
	**/
	int make_gnuplot_charts(std::string const & baseFName,
	                        MatrixD const * const p2forecast,
                            std::string const & gnuplotExe = "gnuplot");
};


/// the main generator
class FcErrTreeGen {
private:
	DimT N; ///< dimension, i.e., the number of variables/margins
	DimT S; ///< number of scenarios
	DimT T; ///< number of periods

	/// \name copula-selection parameters
	/**
		These parameters control which of the possible 2D copulas will be
		used for scenario generation. The more copulas we use, the worse match
		per copula can one expect .. so one should pick the most important ones.

		We provide two methods:
		 - list of pairs
		 - maxDt-variables
		The latter is used if the list is empty

		\todo Add weights for the copulas in the heuristic?
	**/
	///@{
		/// forecasts for the same variable, with varying forecast length
		/**
			For a variable representing (i, dt), we generate 2D-copulas joining
			it with (i, dt ± u), for u = 1..perVarDt.
			perVarDt >= 1
		**/
		int perVarDt;

		/// forecasts for different variables, with varying forecast length
		/**
			for a variable representing (i, dt), we generate 2D-copulas joining
			it with (j, dt ± u), for j≠i and u = 0..intVarDt.
			intVarDt >= 0
		**/
		int intVarDt;

		/// maxDt for all variable pairs
		/**
			If copMaxDt[i, j] = M, then, for a variable representing (i, dt),
			we generate 2D-copulas joining it with (j, dt ± u) for 0 ≤ u ≤ M.

			Note that \c perVarDt corresponds to diagonal values and
			\c intVarDt to off-diagonal values of the matrix
		**/
		ublas::symmetric_matrix<int> copMaxDt;

		/// list of variable pairs in the copula
		std::list<std::pair<DimT, DimT>> copVarPairs;
	///@}

	/// \name main data objects
	/**
		We want to allow both internal data and reference to the caller's data.
		For this, each matrix X is defined as \c MatrixD _X and \c MatrixD & X.
		If we get reference to the caller's data, we simply point \c X to it.
		If we use internal data, these are stored in \c _X and \c X = \c _X.
		In any case, we can simply use \c X throughout the class, as it
		always points to the target data.

		\note The references are denoted as const, so we can use const values
		      in the constructor.

		\warning The above logic works only if the decision is made in the
		         constructor, because C++ references cannot be rebound.
		         However, we want to be able to set the forecast at a later
		         stage -> we use the C++11 reference_wrapper instead
		         of a reference. For now we do this only for the forecast...
	**/
	///@{
		//bool errIsRel;     ///< if true, the errors are relative (in per cent)
		MatrixD _histFErr;  ///< internal version of historical forecast errors; [N, nPts]
		MatrixD const & histFErr = _histFErr; ///< ref. to historical forecast errors
		FcErrTypeList fErrTypes;  ///< forecast-error type of each variable; [N]

		MatrixD _curFcast;  ///< internal version of the current forecast, [T, N]
		/// reference to current forecast, [T, N]
		std::reference_wrapper<const MatrixD> curFcast = std::cref(_curFcast);
	///@}

	bool zeroErrorMeans;  /// < should we fix the error means to zero?

	/// add one stage to an existing scenario tree
	/**
		\param[in]         nBr  number of branches to create at each scenario
		\param[in]          dT  number of periods to add to the tree (>=1)
		\param[out]   copRanks  the generated copula ranks, [nVar, nSc]
		\param[out]  totErrors  converted output values, [T * N, nSc]
		\param[in] p2prevRanks  copula ranks from the prev. iteration, [nVar, nSc]
		//\param[in]  p2forecast  vector of forecasts for the whole horizon; [T, N]
		//\param[out]  p2outTree  the resulting scenario tree (if required)

		// the output scenario tree is needed only at the end
		// in previous iterations, the last parameters would be empty
	**/
	void add_one_stage(DimT const nBr, DimT const dT,
	                   MatrixI const & prevRanks, MatrixI & copRanks,
	                   MatrixD & totErrors);

	/// filter copVars to include only variables up to tMax
	std::list<std::pair<DimT, DimT>> cop_vars_sublist(DimT tMax);

	/// minimal constructor, used only for delegating
	/**
		\param[in]     nmbVars  number of stoch. variables
		\param[in]    histData  historical values and forecasts [nPts, nVars];
		\param[in]  dataFormat  sorting/structure of \c histData
		\param[in]    errTypes  list of error types [nVars]
		\param[in]   zeroMeans  enforce zero means for all errors?
	**/
	FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
	             HistDataFormat const dataFormat,
	             FcErrTypeList const & errTypes,
	             bool const zeroMeans = false);

	/// optional pointer to post-processing of error scenarios (callback)
	ScenPostProcess const * errScenPP = nullptr;

public:
	FcErrTreeGen() = delete; ///< no default/empty constructor

	/// constructor one error type (same for all variables) and max*VarDt
	/**
		\param[in]     nmbVars  number of stoch. variables
		\param[in]    histData  historical values and forecasts [nPts, nVars];
		\param[in]    dataSort  sorting/structure of \c histData
		\param[in]     errType  error type, used for all variables
		\param[in] maxPerVarDt  for var. (i, dt), we generate 2D-copulas with
		                        (i, dt ± u) for u = 1..perVarDt
		\param[in] maxIntVarDt  for var. (i, dt), we generate 2D-copulas with
		                        (j, dt ± u) for u = 0..intVarDt
		\param[in]   zeroMeans  enforce zero means for all errors?
	**/
	FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
	             HistDataFormat const dataFormat,
	             FcErrorType const errType = FcErrorType::linDiff,
	             int const maxPerVarDt = 1, int const maxIntVarDt = 0,
	             bool const zeroMeans = false);

	/// constructor variable error types and specified copula pairs
	/**
		\param[in]     nmbVars  number of stoch. variables
		\param[in]    histData  historical values and forecasts [nPts, nVars];
		\param[in]    dataSort  sorting/structure of \c histData
		\param[in]    errTypes  list of error types [nVars]
		\param[in]      maxDts  matrix of max-dt for all variable pairs
		\param[in]   zeroMeans  enforce zero means for all errors?
	**/
	FcErrTreeGen(DimT const nmbVars, MatrixD const & histData,
	             HistDataFormat const dataFormat,
	             FcErrTypeList const & errTypes,
	             ublas::symmetric_matrix<int> const & maxDts,
	             bool const zeroMeans = false);

	/// set a forecast - this makes a copy of the matrix
	/**
		\param[in] forecast  the current forecast
	**/
	void set_forecast(MatrixD const & forecast) { _curFcast = forecast; }

	/// set a forecast - this uses a reference to the provided matrix
	/**
		\param[in] forecast  the current forecast
	**/
	void set_fcastRef(MatrixD const & forecast) { curFcast = cref(forecast); }

	/// set the list of variable pairs for 2D-copulas
	/**
		This has precedence of MaxDt-type settings

		NB: This operates on variables expended by periods,
		    so it is only valid for a given number of periods!

		\param[in] copVars  list of variable pairs for 2D-copulas
	**/
	void set_cop_var_pairs(std::list<std::pair<DimT, DimT>> const & copVars)
	{
		copVarPairs = copVars;
	}

	/// generate a 2-stage tree (a fan), with given number of scenarios
	/**
		\param[in]    nScen  number of scenarios to generate
		\param[out] outTree  the resulting scenario tree
	**/
	void gen_2stage_tree(DimT const nScen, ScenTree & outTree);

	/// generate a regular multistage tree, with given branching per period
	/**
		\param[in] branching  branching factor per period
		\param[out]  outTree  the resulting scenario tree

		\note By a regular tree we mean a scenario tree where all nodes
		      in the same period have the same number of child nodes.
	**/
	void gen_reg_tree(VectorI const & branching, ScenTree & outTree);

	/// attach a post-processor for error scenarios
	void attach_err_scen_postprocess(ScenPostProcess const * errScPP)
	{ errScenPP = errScPP; }
};

} // namespace FcErr_Gen


// generate logical operators for this enum class
// uses bitmask_operators.hpp, see
// https://www.justsoftwaresolutions.co.uk/cplusplus/using-enum-classes-as-bitfields.html
// NB: must be in the main namespace! for an alternative approach, see
//     http://stackoverflow.com/a/28423827/842693
template<>
struct enable_bitmask_operators<FcErr_Gen::HistDataFormat>{
    static const bool enable=true;
};


#endif // COP_GEN_FC_ERR_HPP
