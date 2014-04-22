#ifndef COP_2D_INFO_HPP
#define COP_2D_INFO_HPP

#include "bivariate_student.h"

// QuantLib libraries used for the normal copula
#include <ql/math/distributions/bivariatenormaldistribution.hpp>
#include <ql/math/distributions/studenttdistribution.hpp>

#include "common.hpp"


// forward declarations, needed by the following class
namespace CopulaDef {
	class CopInfoBy2D;
	class CopInfoData;
	MatrixD & cop_info_data_vals(CopInfoData & copInfo);
	MatrixI & cop_info_data_ranks(CopInfoData & copInfo);
}

namespace Copula2D{

/// specifications of a bivariate copula, used as targets
class Cop2DInfo {
	friend class Cop2DInfTr;
private:

protected:
	/// calculation of one cdf; used mainly to fill the grid
	virtual double calc_cdf(double const u, double const v) const = 0;

	/// \name cdf-grid-related members
	/**
		For some copulas, the calculation of cdf(u,v) might be quite slow.
		For this reason, the main algorithm pre-computes the values and stores
		them in a matrix, to avoid repeating calls. But this will still mean
		repeating calculation in cases where we use the same Cop2DInfo object
		for several pairs of variables. For this reason, it is better to store
		the values inside the class...
	**/
	///@{
		bool useGrid;       ///< use cdf grid for calculations?
		bool customGridPts; ///< custom/non-standard grid positions
		DimT gridN;         ///< size of the grid
		VectorD gridPts;    ///< points of the grid - values from [0,1]
		MatrixF gridCdf;    ///< the computed cdf values (values from [0,1])

		/// convert a value from [0,1] into the position on the grid
		/// \warning at the moment, does not work if customGridPts = true!
		DimT u_to_grid(double const u) const;

		/// this fills the grid with cdf values
		/** assumes values in gridPts; allocates gridCdf if needed **/
		virtual void calc_all_grid_cdfs();
	///@}

public:
	Cop2DInfo()
	: useGrid(false), customGridPts(false), gridN(0) {}
	//{ std::cout << "inside Cop2DInfo()"<< std::endl; }

	virtual ~Cop2DInfo() {}
	//{ std::cout << "inside ~Cop2DInfo(), with gridCdf " << gridCdf.size1() << "x" << gridCdf.size2() << std::endl; }

	/// smart pointer to the class
	typedef boost::shared_ptr<Cop2DInfo> Ptr;

	/// \name cdf-grid- related methods
	///@{
		/// cdf - uses the grid if available, calc_cdf(u,v) otherwise
		virtual double cdf(double const u, double const v) const;

		/// grid cdf, directly using indices of the grid (fastest)
		/**
			\warning no checking -> will fail on index errors!
		**/
		double grid_cdf(DimT const i, DimT const j) const {
			return gridCdf(i, j);
		}

		/// initialize the grid cdf, regular intervals
		/**
			\param[in] N size of the grid
			\param[in] posInInt position of the grid points inside each interval
		**/
		virtual void init_cdf_grid(DimT const N, double const posInInt = 1.0);

		/// initialize the grid cdf, custom grid points
		/**
			\param[in] gridPos  vector of position of the grid points
			\note Since we use custom grid points, ::u_to_grid won't work and
			      hence ::cdf fails -> have to use ::grid_cdf instead!
		**/
		virtual void init_cdf_grid(VectorD const & gridPos);

		MatrixF const & get_cdf_grid() const { return gridCdf; }

		/// delete the grid, to free some memory
		/**
			\todo improve the implementation, add a usage counter

			\warning This is a temporary fix, to decrease memory usage.
				It will create problems if the object is shared by several margin
				pairs - though this is currently implemented only for indep. cop.
				To avoid it, we should add a usage counter - or maybe use the
				counter provided by boost::shared_ptr?
		**/
		void clear_cdf_grid();
	///@}

	//void attach_multivar_info(CopulaDef::CopInfoBy2D * const p2tg,
	//                          DimT const i = -1, DimT const j = -1);
};


/// class for the transposed versions of all copulas
/**
	Logically, this should should not be a derived class, but a simple class
	with a pointer and the two cdf methods. We use inheritance so we can pass
	a pointer as Cop2DInfo*.

	\todo try changing to private or protected inheritance
**/
class Cop2DInfTr : public Cop2DInfo {
private:
	Cop2DInfo const *p2copInfo; ///< non-mutable pointer to the transposed obj.

protected:
	inline double calc_cdf(double const u, double const v) const {
		return p2copInfo->calc_cdf(u, v);
	}
public:
	/// constructor using a (raw) pointer to the transposed object
	Cop2DInfTr(Cop2DInfo const & p2Transp)
	: Cop2DInfo(), p2copInfo(&p2Transp) {}

	/// constructor using a smart pointer to the transposed object
	// the .get() method returns a raw pointer from a smart one
	Cop2DInfTr(Cop2DInfo::Ptr const p2Transp)
	: Cop2DInfo(), p2copInfo(p2Transp.get()) {}
};


//-----------------------------------------------------------------------
//  DERIVED CLASSES

/// The independent copula
class Cop2DIndep : public Cop2DInfo {
	double calc_cdf(const double u, const double v) const { return u * v; }

public:
	Cop2DIndep() : Cop2DInfo() {}
};


/// The Clayton copula
/// copula 4.2.1 from Nelsen, pp. 116
class Cop2DClayton : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf)

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DClayton(const double theta);
};


/// copula 4.2.2 from Nelsen, pp. 116
class Cop2DNelsen2 : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [1, inf)

	double calc_cdf(const double u, const double v) const {
		return std::max(1 - pow(pow(1 - u, th) + pow(1 - v, th), 1.0 / th), 0.0);
	}

public:
	Cop2DNelsen2(const double theta);
};


/// copula 4.2.18 from Nelsen, pp. 118
class Cop2DNelsen18 : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf)

	double calc_cdf(const double u, const double v) const {
		return std::max(1 + th / log(exp(th / (u - 1)) + exp(th / (v - 1))), 0.0);
	}

public:
	Cop2DNelsen18(const double theta);
};


/// Marshall-Olkin copula from Nelsen, pp. 53
/**
	Example of a non-exchangeable copula with singular component.
	The book includes scatter plots for cases (alpha, beta) = (1/2, 3/4) and
	(alpha, beta) = (1/3, 1/4).

	Special cases:
	- alpha = beta -> Caudras-Augé family
	- alpha = 0 or beta = 0 -> product (independence) copula
	- alpha = beta = 1 -> Fréchet-Hoeffding upper bound copula M(u,v) = min(u,v)
**/
class Cop2DMarshallOlkin : public Cop2DInfo {
private:
	double alpha; ///< parameter alpha of the copula; alpha in [0, 1]
	double beta;  ///< parameter beta of the copula; beta in [0, 1]

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DMarshallOlkin(double const alpha_, double const beta_);
};


/// Copula given by a 2D sample
class Cop2DData : public Cop2DInfo {
private:
	/// poiner to the multivar info object
	CopulaDef::CopInfoData * p2multivarTg;
	DimT marg1idx; ///< index of the first margin in the multivar info object
	DimT marg2idx; ///< index of the first margin in the multivar info object

	// have to have it, as it is pure virtual in the base - only a dummy
	double calc_cdf(double const u, double const v) const {
		throw std::logic_error("calc_cdf() should never be used in this class");
		return 0;
	}

	void calc_all_grid_cdfs(); ///< overwriting the base method

protected:
	// - using matrix_row, as this does not allocate new data!
	// - matrix_row does not have a default (empty) constructor, so we cannot
	//   have an array - how would we initialize it?
	ublas::matrix_row<MatrixD> margin1; ///< values of the first margin
	ublas::matrix_row<MatrixD> margin2; ///< values of the first margin
	DimT nPts;     ///< number of the sample/data points in each margin


public:
	/// constructor with a matrix and row numbers
	Cop2DData(MatrixD & histData, int const i, int const j,
	          CopulaDef::CopInfoData * const  multiCopInf = NULL);

	// at the moment, this is not supported
	// would probably have to add a private member with a matrix, since
	// matrix_row need a matrix as a parameter for the constructor
	//Cop2DData(VectorD const & marg1, VectorD const & marg2, int const nSamplPts,
	//          int const gridSize = 0);
};


/// Normal copula
/**
	\note I have not found any library providing the copula cdf, so I use the
	      bivariate distribution and then inverse of the marginal cdfs.
	      For the former, the only implementation with suitable license is from
	      QuantLib; for the latter, we use QuantLib, but boost have the function
	      as well...
**/
class Cop2DNormal : public Cop2DInfo {
private:
	double correl;  ///< correlation

	/// \name QuantLib objects that compute normal CDFs etc.
	///@{
		/// inverse cdf (by default standard normal distrib.)
		QuantLib::InverseCumulativeNormal N01InvCdf;
		/// bivariate cdf (standardized, needs correlation)
		QuantLib::BivariateCumulativeNormalDistribution N01Cdf2D;
	///@}

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DNormal(double const rho);
};


/// Bivariate Student-t copula
/**
	I have not found any library (with suitable license) providing bivariate
	student copula, nor bivariate student distribution. As a result, I had
	implemented the distribution myself and then compute the copula cdf()
	using the univariate inverse CDF from QuantLib.
	(There is also one in boost, but we are already using QuantLib for other
	things...)
**/
class Cop2DStudent : public Cop2DInfo {
private:
	double correl;  ///< correlation
	unsigned dof;   ///< degree of freedom

	///class providing bivariate student cdf
	BivariateCumulativeStudentDistribution tCdf2D;

	/// QuantLib object for the inverse cdf
	QuantLib::InverseCumulativeStudent tInvCdf;

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DStudent(unsigned degF, double const rho);
};

} // namespace


#endif
