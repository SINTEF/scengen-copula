#ifndef COP_SAMPLE_HPP
#define COP_SAMPLE_HPP

//#include <string>

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
private:
	std::vector<bool> haveSc4Marg;

	/// pointer to the object with the multivariate copula info
	/**
		This is used only to adjust probabilities at the very end.
	**/
	CopulaDef::CopulaInfo const * p2copInfo;

	/// array of pointers to the Cop2dInfo objects
	/**
		all (i,j) members with j<i should be of the Cop2DInfTr type,
		pointing to the (j,i)
	**/
	boost::multi_array<Cop2DInfo const*, 2> p2tgInfo;

	/// array of pointers to the Cop2DSample objects
	boost::multi_array<Cop2DSample *, 2> p2sample2D;

	double const *p2prob; ///< pointer to scen. probabilities (can be NULL)

	/// This is the main matrix, including all the scenarios.
	/// Why don't we use IMatrix - is it because we need to get the columns?
	std::vector<IVector> sample;

	/// used to control the level of stochasticity of the results
	double minNumCandPtsRel;  // minimal number of candidate points (in % of all)

protected:
	int nVar; ///< number of random variables
	int nSc;  ///< number of random variables

	double gen_new_margin(int const iNew);

public:
	CopulaSample(int const dim, int const S, double const numCandPtsRel = 0.0);

	virtual ~CopulaSample(void) {};

	/// attach one target 2D copula
	void attach_tg_2Dcop(Cop2DInfo const* p2cop, int const i, int const j,
	                     bool makeTranspTg = true);

	/// attach an external vector of scenario probabilities (no mem. allocation)
	void attach_sc_prob(double const * const scProbs) { p2prob = scProbs; }

	/// attach the object with multivariate copula info
	void attach_tg_cop_info(CopulaDef::CopulaInfo const * const copInf) {
		p2copInfo = copInf;
	}

	/// the main routine; returns the KS-distance
	virtual double gen_sample();

	/// print one 2D copula in a text format
	void print_2D_as_txt(int const i, int const j,
	                     string const fName, bool const sortByScen = false)
	{ p2sample2D[i][j]->print_as_txt(fName, sortByScen); }

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

	int tmp_get_res(int const i, int const s) { return sample[i][s]; }

};

} // namespace

#endif

