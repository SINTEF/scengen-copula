#ifndef COPULA_INFO_HPP
#define COPULA_INFO_HPP

#include <boost/numeric/ublas/symmetric.hpp>  // for correlation matrices
#include <boost/numeric/ublas/triangular.hpp> // pts. to target 2D copulas
#include <map>

#include "common.hpp"
#include "cop2Dinfo.hpp"


namespace CopulaDef{

/// \name objects for the copula name map, used in the main code
///@{
	/// enum for the known copula types
	enum class CopTypeID {sample, normal, indep, student, mixed, unknown};

	/// type for the copula map
	typedef std::map<std::string, CopTypeID> CopNameMapT;

	/// this fills the copula map with the known copula types
	void make_cop_name_map(CopNameMapT & cMap);
///@}


/// class describing the target copula (multivariate)
class CopulaInfo {
private:

protected:
	DimT nVars;   ///< number of random variables
	bool hasCdf;  ///< do we have a formula for multivariate cdf?

	/// \name misc methods used by derived classes
	///@{
		void get_correl_matrix_from_stream(std::istream & is,
		                                   ublas::symmetric_matrix<double> & X);
	///@}

public:
	CopulaInfo(DimT const N, bool const hasMvarCdf)
	: nVars(N), hasCdf(hasMvarCdf) {}

	virtual ~CopulaInfo() {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	/** should throw an error if used on classes with ::haveCdf == false **/
	virtual double cdf(VectorD const u) const = 0;

	DimT dim() const { return nVars; }      ///< get dimension of the copula
	bool has_cdf() const { return hasCdf; } ///< check if we have ::cdf()

	typedef boost::shared_ptr<CopulaInfo> Ptr; ///< smart pointer to the class
};


/*
// ----------------------------------------------------------------------------
/// class for multivariate independent copula
class CopIndep : public CopulaInfo {
public:
	CopIndep(DimT const N)
	: CopulaInfo(N, true) {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	double cdf(VectorD const u) const;

	DimT dim() const { return nVars; }      ///< get dimension of the copula
	bool has_cdf() const { return hasCdf; } ///< check if we have ::cdf()

	typedef boost::shared_ptr<CopulaInfo> Ptr; ///< smart pointer to the class
};
*/


// ----------------------------------------------------------------------------
/// class for copulas specified pairwise by their 2D copulas
class CopInfoBy2D : public CopulaInfo {
public:
	/// type for matrix/array of pointers to the 2D copula info objects
	/**
		When adding a new margin, the generation code connects the new margin
		\c j with the already-generated margins \c i < j, using 2D copulas
		\c (i,j). If follows that we only need a (strict) upper-triangular
		matrix to store the 2D copulas.
	**/
	typedef ublas::triangular_matrix<Copula2D::Cop2DInfo::Ptr,
	                                 ublas::strict_upper> Cop2DInfoPtMatrix;
protected:
	int n2Dcops; ///< number of bivariate copulas

	/// matrix of pointers to the 2D copula specifications
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
		n2Dcops++;
	}

	/// attach one 2D copula info
	/**
		Target given as a shared pointer -> use this for targets that
		already have a shared pointer (i.e. exist somewhere else)
	**/
	void attach_2d_target(Copula2D::Cop2DInfo::Ptr p2tg,
	                      DimT const i, DimT const j) {
		p2Info2D(i, j) = p2tg;
		n2Dcops++;
	}

	Cop2DInfoPtMatrix & get_pts_to_2d_targets() {
		return p2Info2D;
	}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(VectorD const u) const = 0;


	/// initialize cdf grids for all the target 2D copulas; regular intervals
	/**
		\param[in] N size of the grid
		\param[in] useTgPos let the target copula decide the value of \a posInInt
		\param[in] posInInt position of the grid points inside each interval

	**/
	//void init_cdf_grids(DimT const N, bool const useTgPos = true,
	//                    double const posInInt = 0.5);

	/// initialize cdf grids for all the target 2D copulas; custom grid points
	/**
		\param[in] gridPos  vector of position of the grid points
	**/
	//void init_cdf_grids(VectorD const & gridPos);


	/// smart pointer to the class
	typedef boost::shared_ptr<CopInfoBy2D> Ptr;
};


// ----------------------------------------------------------------------------
/// class for multivariate independent copula
class CopInfoIndep : public CopInfoBy2D {
private:
	Copula2D::Cop2DIndep::Ptr p2bivarIndep;

	/// fill #p2Info2D by pointers to #p2bivarIndep; called from constructors
	void setup_2d_targets();

public:
	/// constructor with given dimension
	CopInfoIndep(DimT const N);

	/// constructor with file name of the target distribution
	CopInfoIndep(std::string const & tgFName);

	/// cdf at vector \a u = (u_1, .. , u_n)
	double cdf(VectorD const u) const;

	typedef boost::shared_ptr<CopulaInfo> Ptr; ///< smart pointer to the class
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
		MatrixD hData;  ///< hist. data, original values
		MatrixI hRanks; ///< hist. data, ranks (values from 1 to N)
		MatrixD hU01;   ///< hist. data, ranks scaled to U(0,1)
	///@}

protected:
	/// read the target distribution from a file
	/// \note remember to enclose this in a try{} block!
	void read_tg_file(std::string const & tgFName);

	/// fill \a hRanks and \a hU01 matrices
	void fill_ranks_etc();

	/// creates objects for the 2D targets; called from the constructors
	void setup_2d_targets();

public:
	/// constructor with the target data as input (matrix is [nVar * nPts])
	CopInfoData(MatrixI const & hDataMat);

	/// constructor with file name of the target distribution
	CopInfoData(std::string const & tgFName);

	/// constructor with the matrix of ranks as input
	//CopInfoData(TMatrixI & ranksMat, bool const fillU01Data = true);

	virtual ~CopInfoData() {}

	MatrixD & data_vals() { return hData; }
	MatrixI & data_ranks() { return hRanks; }
	MatrixD & data_u01() { return hU01; }

	double cdf(VectorD const u) const;
};

// non-member accessors to data of the CopInfoData class - used because we
// cannot forward-declare the members from within cop2Dinfo.hpp!
MatrixD & cop_info_data_vals(CopInfoData & copInfo);
MatrixI & cop_info_data_ranks(CopInfoData & copInfo);


// ----------------------------------------------------------------------------
/// normal copula, i.e. a collection of bivariate normal copulas
class CopInfoNormal : public CopInfoBy2D {
private:

protected:
	ublas::symmetric_matrix<double> correlMat; ///< correlation matrix

	/// read the correlation matrix from a file
	/// \note remember to enclose this in a try{} block!
	void read_correl_mat(std::string const & tgFName);

	/// creates objects for the 2D targets; called from the constructors
	virtual void setup_2d_targets();

	/// default constructor - needed from derived classes
	/*
		This does not work, as it needs an empty constructor of CopInfoBy2D,
		which does not exist.
		A better option could be to have constructors using open input streams
		and external functions to create objects based on a file name.
		This way, we could read some parts of the file before calling the
		constructors, something not possible from within the class(?)
		Notes: - the functions should/could be templated
		       - these constructors could be protected, or perhaps even private,
		         with the external functions declared as friends?
	*/
	CopInfoNormal() : CopInfoBy2D(0) {}

public:
	/// constructor with the target data as input
	CopInfoNormal(MatrixD const & correls);

	/// constructor with an open input stream
	//CopInfoNormal(istream & is);

	/// constructor with file name of the target distribution
	CopInfoNormal(std::string const & tgFName);

	virtual ~CopInfoNormal() {}

	double cdf(VectorD const u) const {
		throw std::logic_error("class CopInfoNormal does not have cdf()");
	}
};


// ----------------------------------------------------------------------------
/// student-t copula, i.e. a collection of bivariate student-t copulas
/**
	\note Deriving from the normal copula, to get some of the methods
**/
class CopInfoStudent : public CopInfoNormal {
private:
	unsigned dof;                              ///< degree of freedom

protected:
	/// creates objects for the 2D targets; called from the constructors
	void setup_2d_targets() override;

	/// read the parameters (dof and correlation matrix) from a file
	/// \note remember to enclose this in a try{} block!
	void get_params_from_file(std::string const & tgFName);

public:
	/// constructor with the target data as input
	CopInfoStudent(unsigned degF, MatrixD const & correls);

	/// constructor with dof and name of file with the correlation matrix
	CopInfoStudent(unsigned degF, std::string const & tgFName);

	/// constructor with name of the parameter file
	/**
		the first number in the file is interpreted as the dof,
		the rest as the correlation matrix
	**/
	CopInfoStudent(std::string const & tgFName);

	virtual ~CopInfoStudent() {}

	double cdf(VectorD const u) const {
		throw std::logic_error("class CopInfoStudent does not have cdf()");
	}
};


// ----------------------------------------------------------------------------
/// target copula given as a list of 2D copulas
class CopInfoGen2D : public CopInfoBy2D {

public:
	/// constructor with the input file name as a param
	CopInfoGen2D(std::string const & tgFName);

	double cdf(VectorD const u) const {
		throw std::logic_error("class CopInfoGen2D does not have cdf()");
	}
};


} // namespace

#endif
