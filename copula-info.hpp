#ifndef COPULA_INFO_HPP
#define COPULA_INFO_HPP

#include <boost/numeric/ublas/triangular.hpp>

#include "common.hpp"
#include "cop2Dinfo.hpp"


namespace CopulaDef{

/// class describing the target copula (multivariate)
class CopulaInfo {
private:

protected:
	DimT nVars;   ///< number of random variables
	bool hasCdf; ///< do we have a formula for multivariate cdf?

public:
	CopulaInfo(DimT const N, bool const hasMvarCdf)
	: nVars(N), hasCdf(hasMvarCdf) {}

	virtual ~CopulaInfo() {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	/** should throw an error if used on classes with ::haveCdf == false **/
	virtual double cdf(UVector const u) const = 0;

	DimT dim() const { return nVars; }      ///< get dimension of the copula
	bool has_cdf() const { return hasCdf; } ///< check if we have ::cdf()

	typedef boost::shared_ptr<CopulaInfo> Ptr; ///< smart pointer to the class
};


// ----------------------------------------------------------------------------
/// class for copulas specified pairwise by their 2D copulas
class CopInfoBy2D : public CopulaInfo {
public:
	typedef ublas::triangular_matrix<Copula2D::Cop2DInfo::Ptr,
	                                 ublas::strict_upper> Cop2DInfoPtMatrix;
protected:
	int n2Dcops; ///< number of bivariate copulas

	/// matrix of pointers to the 2D copula specifications
	/**
		When adding a new margin, the generation code connects the new margin
		\c j with the already-generated margins \c i < j, using 2D copulas
		\c (i,j). If follows that we only need a (strict) upper-triangular
		matrix to store the 2D copulas.
	**/
	Cop2DInfoPtMatrix p2Info2D;
	//boost::multi_array<Copula2D::Cop2DInfo::Ptr, 2> p2Info2D;

public:
	CopInfoBy2D(DimT const N, bool const hasCdf = false)
	: CopulaInfo(N, hasCdf), n2Dcops(0), p2Info2D()
	{
		// resize p2Info2D - cannot be done in the initialization, since
		// a strict triangular matrix fails if initialized with (0,0)!
		if (N > 0) {
			p2Info2D.resize(N, N);
		}
	}

	virtual ~CopInfoBy2D() {}

	/// attach one 2D copula info
	/**
		Target given as a 'raw' pointer -> use this for targets that are created
		just for calling this function (and hence do not exist anywhere else).
	**/
	void attach_2d_target(Copula2D::Cop2DInfo * p2tg,
	                      DimT const i, DimT const j) {
		p2Info2D(i, j).reset(p2tg);
	}

	/// attach one 2D copula info
	/**
		Target given as a shared pointer -> use this for targets that
		already have a shared pointer (i.e. exist somewhere else)
	**/
	void attach_2d_target(Copula2D::Cop2DInfo::Ptr p2tg,
	                      DimT const i, DimT const j) {
		p2Info2D(i, j) = p2tg;
	}

	Cop2DInfoPtMatrix & get_pts_to_2d_targets() {
		return p2Info2D;
	}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(UVector const u) const = 0;

	typedef boost::shared_ptr<CopInfoBy2D> Ptr;
};


// ----------------------------------------------------------------------------
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
	void read_tg_file(std::string const & tgFName);

	/// fill \a hRanks and \a hU01 matrices
	void fill_ranks_etc();

public:
	/// constructor with the target data as input
	CopInfoData(UIMatrix const & hDataMat);

	/// constructor with file name of the target distribution
	CopInfoData(std::string const & tgFName);

	/// constructor with the matrix of ranks as input
	//CopInfoData(TMatrixI & ranksMat, bool const fillU01Data = true);

	virtual ~CopInfoData() {}

	/// creates objects for the 2D targets
	void setup_2d_targets();

	UMatrix & data_vals() { return hData; }
	UIMatrix & data_ranks() { return hRanks; }
	UMatrix & data_u01() { return hU01; }

	double cdf(UVector const u) const;
};

// non-member accessors to data of the CopInfoData class - used because we
// cannot forward-declare the members from within cop2Dinfo.hpp!
UMatrix & cop_info_data_vals(CopInfoData & copInfo);
UIMatrix & cop_info_data_ranks(CopInfoData & copInfo);


// ----------------------------------------------------------------------------
/// target copula described by a historical data (sample)
class CopInfoNormal : public CopInfoBy2D {
private:
	UMatrix correlMat; ///< correlation matrix

protected:
	/// read the correlation matrix from a file
	/// \note remember to encluse this in a try{} block!
	void read_correl_mat(std::string const & tgFName);

public:
	/// constructor with the target data as input
	CopInfoNormal(UMatrix const & correls);

	/// constructor with file name of the target distribution
	CopInfoNormal(std::string const & tgFName);

	virtual ~CopInfoNormal() {}

	/// creates objects for the 2D targets
	void setup_2d_targets();

	double cdf(UVector const u) const {
		throw std::logic_error("class CopInfoNormal does not have cdf()");
	}
};


} // namespace

#endif
