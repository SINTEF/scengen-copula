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
	DimT nVars;      ///< number of random variables

public:
	CopulaInfo(DimT const N = 0) : nVars(N) {}

	virtual ~CopulaInfo() {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(TVectorD const u) const = 0;

	DimT dim() const { return nVars; }

	typedef boost::shared_ptr<CopulaInfo> Ptr;
};


/// class for copulas specified pairwise by their 2D copulas
class CopInfoBy2D : public CopulaInfo {
protected:
	int n2Dcops; ///< number of bivariate copulas
	///< matrix of pointers to the 2D copula specifications
	boost::multi_array<Copula2D::Cop2DInfo::Ptr, 2> p2Info2D;

public:
	CopInfoBy2D(int const N = 0)
	: CopulaInfo(N), n2Dcops(0), p2Info2D(boost::extents[N][N])
	{}

	virtual ~CopInfoBy2D() {}

	/// attach one 2D copula info
	/**
		Target given as a 'raw' pointer -> use this for targets that are created
		just for calling this function (and hence do not exist anywhere else).
	**/
	void attach_2d_target(Copula2D::Cop2DInfo * p2tg, int const i,
	                      int const j, bool const makeTransp = true) {
		p2Info2D[i][j].reset(p2tg);
		if (makeTransp) {
			assert (p2Info2D[j][i] == false);
			p2Info2D[j][i].reset(new Copula2D::Cop2DInfTr(p2Info2D[i][j]));
		}
	}

	/// attach one 2D copula info
	/**
		Target given as a shared pointer -> use this for targets that
		already have a shared pointer (i.e. exist somewhere else)
	**/
	void attach_2d_target(Copula2D::Cop2DInfo::Ptr p2tg, int const i,
	                      int const j, bool const makeTransp = true) {
		p2Info2D[i][j] = p2tg;
		if (makeTransp) {
			assert (p2Info2D[j][i] == false);
			p2Info2D[j][i].reset(new Copula2D::Cop2DInfTr(p2Info2D[i][j]));
		}
	}

	boost::multi_array<Copula2D::Cop2DInfo::Ptr, 2> & get_pts_to_2d_targets() {
		return p2Info2D;
	}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(TVectorD const u) const = 0;

	typedef boost::shared_ptr<CopInfoBy2D> Ptr;
};


/// target copula described by a historical data (sample)
class CopInfoData : public CopInfoBy2D {
protected:
	DimT nPts; ///< number of the sample/data points (in the hist. data)

private:
	/// \name historical data, different formats
	/**
		All matrices have margins in rows, i.e. i-th margin is (i,*).
	**/
	///@{
		UMatrix hData;   ///< hist. data, original values
		UIMatrix hRanks; ///< hist. data, ranks (values from 1 to N)
		UMatrix  hU01;   ///< hist. data, ranks scaled to U(0,1)
	///@}

protected:
	/// read the target distribution from a file
	/// \note remember to encluse this in a try{} block!
	void read_tg_file(std::string tgFName);

	/// fill \a hRanks and \a hU01 matrices
	void fill_ranks_etc();

public:
	/// constructor with file name of the target distribution
	CopInfoData(std::string & tgFName);

	/// constructor with the target data as input
	CopInfoData(UIMatrix const & hDataMat);

	/// constructor with the matrix of ranks as input
	//CopInfoData(TMatrixI & ranksMat, bool const fillU01Data = true);

	virtual ~CopInfoData() {}

	/// creates objects for the 2D targets
	/**
		\param[in] makeTranspTgs should we create the transposed targets?
		           (if no, only targets with j>i get created)
	**/
	void setup_2d_targets(bool const makeTranspTgs = true);

	UMatrix & data_vals() { return hData; }
	UIMatrix & data_ranks() { return hRanks; }
	UMatrix & data_u01() { return hU01; }

	double cdf(TVectorD const u) const;
};

// non-member accessors to data of the CopInfoData class - used because we
// cannot forward-declare the members from within cop2Dinfo.hpp!
UMatrix & cop_info_data_vals(CopInfoData & copInfo);
UIMatrix & cop_info_data_ranks(CopInfoData & copInfo);

} // namespace

#endif
