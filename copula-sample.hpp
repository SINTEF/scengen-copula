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

	/// array of pointers to the Cop2DSample objects
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

	/// this sets/fixes values for some margins
	/**
        \param[in] X     matrix with the ranks being fixed [m x nSc]
        \param[in] marg  margin indexes to fix [m]; if empty, fix the first m margins
	**/
	void fix_marg_values(MatrixI const & X, VectorI const & marg = VectorI());

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

	/// get the ranks as one big matrix
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

	//int tmp_get_res(int const i, int const s) { return sample[i][s]; }
};

} // namespace

#endif

