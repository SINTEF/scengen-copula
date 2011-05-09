#ifndef COPULA_INFO_HPP
#define COPULA_INFO_HPP

#include <cassert>
#include "common.hpp"
#include "cop2Dinfo.hpp"

namespace CopulaDef{

/// class describing the target copula (multivariate)
class CopulaInfo {
private:

protected:
	int nVars;      ///< number of random variables

public:
	CopulaInfo(int const N = 0) : nVars(N) {}

	virtual ~CopulaInfo() {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(TVectorD const u) const = 0;
};


/// class for copulas specified pairwise by their 2D copulas
class CopInfoBy2D : public CopulaInfo {
protected:
	int n2Dcops; ///< number of bivariate copulas
	///< matrix of pointers to the 2D copula specifications
	boost::multi_array<Copula2D::Cop2DInfo::Ptr, 2> p2Info2D;

public:
	CopInfoBy2D(int const N = 0) : CopulaInfo(N), n2Dcops(0),
	p2Info2D(boost::extents[N][N]) {}

	virtual ~CopInfoBy2D() {}

	TO DO: write implementation
	       change CopInfoData::setup_2d_targets() to use this function...!!!
	attach_2d_target(Copula2D::Cop2DInfo::Ptr p2tg, int const i, int const j,
	                 bool const makeTransp = true);

	boost::multi_array<Copula2D::Cop2DInfo::Ptr, 2> & get_pts_to_2d_targets() {
		return p2Info2D;
	}
};


/// target copula described by a historical data (sample)
class CopInfoData : public CopInfoBy2D {
protected:
	int nPts; ///< number of the sample/data points

private:
	/// \name historical data, different formats
	/**
		All matrices have margins in rows, i.e. i-th margin is (i,*).
	**/
	///@{
	// we use vector of vectors for hData, so we can easily pass rows as
	// pointers or references - ublas gives access to rows and columns,
	// but only to create new object (copy values)?
	std::vector<UVector> hData; ///< hist. data, original values
	UIMatrix hRanks; ///< hist. data, ranks (values from 1 to N)
	UMatrix  hU01;   ///< hist. data, ranks scaled to U(0,1)
	///@}

protected:
	/// read the target distribution from a file
	/// \note remember to encluse this in a try{} block!
	void read_tg_file(std::string tgFName);

	/// fill \a hRanks and \a hU01 matrices
	void fill_ranks_etc();

	/// creates objects for the 2D targets
	/**
		\param[in] makeTranspTgs should we create the transposed targets?
		           (if no, only targets with j>i get created)
	**/
	void setup_2d_targets(bool const makeTranspTgs = true);

public:
	/// constructor with file name of the target distribution
	CopInfoData(std::string & tgFName);

	/// constructor with the target data as input
	CopInfoData(UIMatrix const & hDataMat);

	/// constructor with the matrix of ranks as input
	//CopInfoData(TMatrixI & ranksMat, bool const fillU01Data = true);

	virtual ~CopInfoData() {}

	double cdf(TVectorD const u) const;
};

} // namespace

#endif
