#ifndef COP_SAMPLE_HPP
#define COP_SAMPLE_HPP

//#include <string>
#include <utility> // to get std::pair<>

#include "common.hpp"
#include "cop2Dinfo.hpp"
#include "cop2Dsample.hpp"
#include "copula-info.hpp"

//using std::string;
using Copula2D::Cop2DInfo;
using Copula2D::Cop2DSample;


namespace CopulaScen {

/// the main class for creating a multivariate sample
class CopulaSample {
protected:
	DimT nVar; ///< number of random variables
	DimT nSc;  ///< number of scenarios

private:
	/// indicator if we have generated scenarios for a given margin
	std::vector<bool> haveSc4Marg;

	/// pointer to the object with the multivariate copula info
	/**
		This is used only to adjust probabilities at the very end.
	**/
	CopulaDef::CopulaInfo const * p2copInfo;

	/// matrix of pointers to the Cop2dInfo objects
	/**
		It is a (strict) upper-triangular matrix, i.e. we have (i,j) with j > i.
	**/
	CopulaDef::CopInfoBy2D::Cop2DInfoPtMatrix & p2tgInfo;

	/// array of pointers to the Cop2DSample objects \todo use a smart pointer!?
	Array2D<Cop2DSample *> p2sample2D;

	double const *p2prob; ///< pointer to scen. probabilities (can be NULL)

	/// This is the main matrix, including all the scenarios.
	/**
		Dimension is [nVar x nSc]

		\note We use a vector of vectors, for easier passing of values for one variable
	**/
	std::vector<VectorI> sample;

	/// used to control the level of stochasticity of the results
	unsigned minNumCandScens;  // minimal number of candidate points

	/// number of margins with allowed duplicate values
	/**
		If we start with fixing the first M margins, these may include
		duplicate values (the same rank in several scenarios).
		Normally, this is not allowed. Setting this value to M will allow
		duplicate values in all margins (i,j) with i<M.
	**/
	DimT nmbMargsWithDuplicates = 0;

protected:
	double gen_new_margin(DimT const iNew);

public:
	CopulaSample(CopulaDef::CopInfoBy2D::Ptr p2tg, DimT const S,
	             unsigned const numCandPts = 1);

	virtual ~CopulaSample() {}

	/// attach an external vector of scenario probabilities (no mem. allocation)
	void attach_sc_prob(double const * const scProbs) { p2prob = scProbs; }

	/// attach the object with multivariate copula info
	void attach_tg_cop_info(CopulaDef::CopulaInfo const * const copInf) {
		p2copInfo = copInf;
	}

	/// the main routine; returns the KS-distance
	virtual double gen_sample();

	/// get the ranks as one big matrix, [nVar, nSc]
	void get_result_ranks(MatrixI & ranks);

	/// print one 2D copula in a text format
	void print_2D_as_txt(int const i, int const j,
	                     string const fName, bool const scaleTo01,
	                     bool const sortByScen = false)
	{ p2sample2D(i,j)->print_as_txt(fName, scaleTo01, sortByScen); }

	/// print the whole matrix as text
	/**
		\param[in] fName name of the output file
		\param[in] scaleTo01 scale from 0..N-1 to (0,1)
		\param[in] sortByMarg sort output by values of a margin; <0 means no sort
	**/
	void print_as_txt(string const fName, bool const scaleTo01,
	                  int const sortByMarg = -1);

	/// write a data file for the probability-allocation model
	void write_gmp_data(string const fName = "prob-alloc.dat");

	/// shuffle the resulting matrix (otherwise, the first margin is sorted)
	void shuffle_results();

	/// this sets/fixes values for the first margins
	/**
        \param[in]               X  matrix with the ranks being fixed; [m x nSc]
        \param[in] allowDuplicates  Do we allow multiple scenarios with the same
                                    rank? Normally, this should be \c false,
                                    but we use it for multi-stage trees
                                    generated from forecast errors.

        \note Must be called before \a gen_sample().
        \note We allow only fixing of the first margins. The reason is that
              if we fixed only margin 2, then generation of margin 1 should
              use target 2D copula(2,1) as well .. but these targets with i>j
              are not created at the moment - and target(j,i) is, in general,
              different from target(i,j)!
	**/
	void fix_marg_values(MatrixI const & X, bool const allowDuplicates = false);

	//int tmp_get_res(int const i, int const s) { return sample[i][s]; }
};

} // namespace

#endif

