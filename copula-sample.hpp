#ifndef COP_SAMPLE_HPP
#define COP_SAMPLE_HPP

//#include <string>

#include "common.hpp"
#include "cop2Dinfo.hpp"
#include "cop2Dsample.hpp"

//using std::string;
using Copula2D::Cop2DInfo;
using Copula2D::Cop2DSample;


namespace CopulaScen {

/// the main class for creating a multivariate sample
class CopulaSample {
private:
	std::vector<bool> haveSc4Marg;

	/// array of pointers to the Cop2dInfo class
	/**
		all (i,j) members with j<i should be of the Cop2DInfTr type,
		pointing to the (j,i)
	**/
	boost::multi_array<Cop2DInfo const*, 2> p2tgInfo;

	boost::multi_array<Cop2DSample *, 2> p2sample2D;

	double const *p2prob; ///< pointer to scen. probabilities (can be NULL)

	/// This is the main matrix, including all the scenarios.
	/// Why don't we use IMatrix - is it because we need to get the columns?
	std::vector<IVector> sample;

protected:
	int nVar; ///< number of random variables
	int nSc;  ///< number of random variables

	double gen_new_margin(int const iNew);

public:
	CopulaSample(int const dim, int const S);

	virtual ~CopulaSample(void) {};

	void attach_tg_2Dcop(Cop2DInfo const* p2cop, int const i, int const j,
											 bool makeTranspTg = true);

	/// the main routine; returns the KS-distance
	virtual double gen_sample();

	/// print one 2D copula in a text format
	void print_2D_as_txt(int const i, int const j,
											 string const fName, bool const sortByScen = false)
	{ p2sample2D[i][j]->print_as_txt(fName, sortByScen); }

	/// print the whole matrix as text
	/**
		\param[in] sortByMarg sort output by values of a margin; <0 means no sort
	**/
	void print_as_txt(string const fName, int const sortByMarg = -1);

};

} // namespace

#endif

