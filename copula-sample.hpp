#ifndef COP_SAMPLE_HPP
#define COP_SAMPLE_HPP

//#include <string>

#include "common.hpp"
#include "cop2Dinfo.hpp"

//using std::string;

namespace Copula2D{

/// the main class for creating a multivariate sample
class CopulaSample {
private:
	std::vector<bool> haveSc4Marg;

	/// array of pointers to the Cop2dInfo class
	/**
		all i>j members should be of the Cop2DInfTr, pointing to the (j,i)
	**/
	multi_array<Cop2DInfo const*, 2> p2tgInfo;

	multi_array<Cop2DSample *, 2> p2sample2D;

	double const *p2prob; ///< pointer to scen. probabilities (can be NULL)

	std::vector<IVector> sample;

protected:
	int nVar; ///< number of random variables
	int nSc;  ///< number of random variables

	double gen_new_margin(int const iNew);

public:
	CopulaSample(int const dim);

	virtual ~CopulaSample() {};

	// ! check that p2cop->N == nS !!!
	void attach_tg_2Dcop(Cop2DInfo const* p2cop, int const i, int const j
											 bool makeTranspTg = true);

	/// the main routine; returns the KS-distance
	virtual double gen_sample();
};
