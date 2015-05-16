#ifndef COP_2D_INFO_HPP
#define COP_2D_INFO_HPP

#include "bivariate_student.h"
#include "common.hpp"

// QuantLib libraries used for the normal copula
#include "external/QuantLib/math/distributions/bivariatenormaldistribution.hpp"
#include "external/QuantLib/math/distributions/bivariatestudenttdistribution.hpp"
#include "external/QuantLib/math/distributions/studenttdistribution.hpp"

#include <boost/optional.hpp>

#include <iostream>
#include <map>


// forward declarations, needed by the following class
namespace CopulaDef {
	class CopInfoBy2D;
	class CopInfoData;
	MatrixD & cop_info_data_vals(CopInfoData & copInfo);
	MatrixI & cop_info_data_ranks(CopInfoData & copInfo);
}

namespace Copula2D{

/// \name objects for the copula name map, used in the main code
///@{
/*
	/// enum for the known copula types
	enum class Cop2DTypeID {indep, Clayton, Gumbel, Frank, Nelsen2, Nelsen18,
	                        MarshallOlkin, data, normal, student};

	/// type for the copula map
	typedef std::map<std::string, Cop2DTypeID> Cop2DNameMapT_Old;

	/// this fills the copula map with the known copula types
	void make_2d_cop_name_map(Cop2DNameMapT_Old & cMap);
*/

	class Cop2DInfo; ///< forward declaration needed for UnivMargNameMapT

	/// map from margin's name to its constructor
	/// \todo return a smart pointer instead??
	typedef std::map<std::string, Cop2DInfo * (*) (std::istream &)>
	Cop2DNameMapT;
///@}


/// the base class of all specifications of bivariate copulas
class Cop2DInfo {
	//friend class Cop2DInfTr;
private:
	/// template for adding new entries to a shape-name-map
	/// \todo checks whether the string is there already
	template <class T>
	static void add_to_map(std::string const & idStr);

	/// the name map
	/**
		converts from margin name to a constructor that takes a string (params
		from the input file) and integer (number of scenarios)

		/// \todo Add all the name-map related objects and methods to a separate
		///       class and derive all classes that use it from it?
	**/
	static Cop2DNameMapT nameMap;

	/// this fills the map
	static void init_name_map();

	/// \name members and methods for converting between ranks and values
	///@{
		void gen_cop_pts(); ///< this generates values for \c copPts

protected:
		DimT nSc;           ///< number of scenarios/samples; needed for \c cdfR()
		VectorD copPts;     ///< vector of points of copula discretization
	///@}


public:
	Cop2DInfo(DimT const nScens = 0)
	: nSc(nScens), copPts(nScens)
	{
		if (nScens > 0)
			gen_cop_pts();
	}

	virtual ~Cop2DInfo() {}

	/// smart pointer to the class
	typedef boost::shared_ptr<Cop2DInfo> Ptr;

	/// getter for the number of scenarios/samples
	DimT nmb_sc() const { return nSc; }

	/// cdf, using percentiles - not used anywhere!
	virtual double cdf(double const u, double const v) const = 0;

	/// cdf using ranks- this is used by the copula-generating alg.
	virtual double cdfR(DimT const i, DimT const v) const = 0;

	/// creates a new object based on its name
	static Cop2DInfo * make_new(std::string const & margName,
	                            std::istream & paramStr);

	virtual void set_nmb_scens(DimT const nScens);
};


/// base class for copula classes that actually perform the calculations
class Cop2DComp : public Cop2DInfo {
private:

protected:
	/// calculation of one cdf; this is what the derived classes must override
	virtual double calc_cdf(double const u, double const v) const = 0;

public:
	Cop2DComp(DimT const nmbScen = 0)
	: Cop2DInfo(nmbScen) {}

	virtual ~Cop2DComp() {}

	/// cdf, using percentiles - not used anywhere!
	virtual double cdf(double const u, double const v) const;

	/// cdf using ranks- this is used by the copula-generating alg.
	virtual double cdfR(DimT const i, DimT const v) const;
};


//-----------------------------------------------------------------------
//  VIEW CLASSES
/// \name Classes providing different 'views' on copulas
///@{

	/// class for the transposed versions of all copulas
	/**
		\todo try changing to private or protected inheritance
	**/
	class Cop2DInfTr : public Cop2DInfo {
	private:
		Cop2DInfo const *p2copInfo; ///< non-mutable pointer to the transposed obj.

	protected:

	public:
		/// constructor using a reference to the transposed object
		Cop2DInfTr(Cop2DInfo const & transpCop)
		: Cop2DInfo(transpCop.nmb_sc()), p2copInfo(&transpCop) {}

		/// constructor using a smart pointer to the transposed object
		// the .get() method returns a raw pointer from a smart one
		Cop2DInfTr(Cop2DInfo::Ptr const p2transp)
		: Cop2DInfo(p2transp->nmb_sc()), p2copInfo(p2transp.get()) {}
	};


	/// copula with all the values in a pre-computed grid
	/**
		This won't save time if we use the object only once, as the algorithm
		already ensures that each cdf is calculated only one.
		In other words, this makes sense only if we have more margins with
		the same copula.

		Storing grid for N scenarios needs N*N double values, so it becomes
		problematic for large N (ca. thousand and more).
		This class addresses this by having a grid with some given size, smaller
		than the number of scenarios, and then computed CDFs by approximation
		from the grid.
		The grid interpolation uses formulas from
		https://en.wikipedia.org/wiki/Bilinear_interpolation.

		\note Unlike the exact grid, the interpolation grid needs values
		      at both extremes, so we need (N + 1) values for N intervals.
	**/
	class Cop2DGrid : public Cop2DInfo {
	private:
		/// \name cdf-grid-related members
		///@{
			/// convert a value from [0,1] into the rank/scenario index
			boost::optional<DimT> perc_to_rank(double const u) const;
		///@}

	protected:
		Cop2DInfo::Ptr const p2copInfo; ///< pointer to the copula we are viewing

		/// \name cdf-grid-related members - size of the grid = nSc
		///@{
			MatrixD gridCdf;    ///< the computed cdf values (values from [0,1])

			// assumes values in gridPts; allocates gridCdf if needed
			virtual void calc_all_grid_cdfs();
		///@}

		/// \name members and methods for interpolated grid
		///@{
			bool exactGrid;       ///< does the grid coincide witch scenarios?
			DimT gridN;           ///< number of grid points
			VectorD gridPts;      ///< grid values
			void gen_grid_pts();  ///< generate the grid values
			void set_grid_size(); ///< find grid size (if not given)

			double cdfMult;      ///< scaling constant used in the interpolation
			double gridMult;     ///< difference between the grid sizes
		///@}


	public:
		/// constructor using a smart pointer to the actual copula
		Cop2DGrid(Cop2DInfo::Ptr const p2tgCop, DimT const gridSize = 0,
		          bool const fillGrid = true)
		: Cop2DInfo((p2tgCop ? p2tgCop->nmb_sc() : 0)),
		  p2copInfo(p2tgCop),
		  exactGrid(gridSize == 0),
		  gridN((exactGrid ? nSc : gridSize)),
		  gridPts((exactGrid ? 0 : gridSize + 1)),
		  cdfMult(pow(gridSize, 2)), gridMult((double) gridSize / nSc)
		{
			set_grid_size();
			if (gridN > 0) {
				if (!exactGrid) {
					gen_grid_pts();
				}
				calc_all_grid_cdfs();
			}
		}

		virtual double cdf(double const u, double const v) const;

		virtual double cdfR(DimT const i, DimT const j) const;

		virtual void set_nmb_scens(DimT const nScens) override
		{
			Cop2DInfo::set_nmb_scens(nScens); // parent-class version
			calc_all_grid_cdfs();
		}
	};


//-----------------------------------------------------------------------
//  DERIVED CLASSES

/// The independent copula
class Cop2DIndep : public Cop2DComp {
	double calc_cdf(const double u, const double v) const { return u * v; }

public:
	Cop2DIndep() : Cop2DComp() {}

	Cop2DIndep(std::istream & paramStr) : Cop2DComp() {}
};


/// The Clayton copula
/// copula 4.2.1 from Nelsen, pp. 116
class Cop2DClayton : public Cop2DComp {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf) - {0}

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DClayton(const double theta);

	Cop2DClayton(std::istream & paramStr);
};


/// The Gumbel copula
/// copula 4.2.4 from Nelsen, pp. 116
class Cop2DGumbel : public Cop2DComp {
private:
	double th;  ///< parameter theta of the copula; th in [1, inf)
	double iTh; ///< 1 / th

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DGumbel(const double theta);

	Cop2DGumbel(std::istream & paramStr);
};


/// The Frank copula
/// copula 4.2.5 from Nelsen, pp. 116
class Cop2DFrank : public Cop2DComp {
private:
	double th; ///< parameter theta of the copula; th in R - {0}
	double C;  ///< exp(-th)

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DFrank(const double theta);

	Cop2DFrank(std::istream & paramStr);
};


/// copula 4.2.2 from Nelsen, pp. 116
class Cop2DNelsen2 : public Cop2DComp {
private:
	double th; ///< parameter theta of the copula; th in [1, inf)

	double calc_cdf(const double u, const double v) const {
		return std::max(1 - pow(pow(1 - u, th) + pow(1 - v, th), 1.0 / th), 0.0);
	}

public:
	Cop2DNelsen2(const double theta);

	Cop2DNelsen2(std::istream & paramStr);
};


/// copula 4.2.18 from Nelsen, pp. 118
class Cop2DNelsen18 : public Cop2DComp {
private:
	double th; ///< parameter theta of the copula; th in [2, inf)

	double calc_cdf(const double u, const double v) const {
		return std::max(1 + th / log(exp(th / (u - 1)) + exp(th / (v - 1))), 0.0);
	}

public:
	Cop2DNelsen18(const double theta);

	Cop2DNelsen18(std::istream & paramStr);
};


/// copula 4.2.21 from Nelsen, pp. 118
class Cop2DNelsen21 : public Cop2DComp {
private:
	double th;  ///< parameter theta of the copula; th in [1, inf)
	double iTh; ///< 1 / th

	double calc_cdf(const double u, const double v) const {
		return 1 - pow(1 - pow(std::max(
				pow(1 - pow(1 - u, th), iTh) + pow(1 - pow(1 - v, th), iTh) - 1,
			0.0), th), iTh);
	}

public:
	Cop2DNelsen21(const double theta);

	Cop2DNelsen21(std::istream & paramStr);
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
class Cop2DMarshallOlkin : public Cop2DComp {
private:
	double alpha; ///< parameter alpha of the copula; alpha in [0, 1]
	double beta;  ///< parameter beta of the copula; beta in [0, 1]

	double calc_cdf(const double u, const double v) const;

public:
	Cop2DMarshallOlkin(double const alpha_, double const beta_);

	Cop2DMarshallOlkin(std::istream & paramStr);
};


/// Copula given by a 2D sample
/**
	\warning
	For the moment, this works only if all copulas are given by historical
	data, since it simply points to the multivariate \c CopInfoData object.
**/
class Cop2DData : public Cop2DGrid {
private:
	/// poiner to the multivar info object
	CopulaDef::CopInfoData * p2multivarTg;
	DimT marg1idx;   ///< index of the first margin in the multivar info object
	DimT marg2idx;   ///< index of the first margin in the multivar info object

	void calc_all_grid_cdfs() override; ///< overwriting the base method

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
	          CopulaDef::CopInfoData * const  p2CopInf = nullptr,
	          DimT nScens = 0);

	// at the moment, this is not supported
	// would probably have to add a private member with a matrix, since
	// matrix_row need a matrix as a parameter for the constructor
	//Cop2DData(VectorD const & marg1, VectorD const & marg2, int const nSamplPts,
	//          int const gridSize = 0);

	virtual void set_nmb_scens(DimT const nScens) override;
};


/// Copula given by a 2D sample, variant for many scenarios
/**
	\warning
	For the moment, this works only if all copulas are given by historical
	data, since it simply points to the multivariate \c CopInfoData object.
**/
/*
class Cop2DBigData : public Cop2DData {
private:
	void adjust_cdf_grid(); ///< add zero rows and columns (for interpolation)

	void calc_all_grid_cdfs() override; ///< overwriting the base method

protected:

public:
	/// constructor with a matrix and row numbers
	Cop2DBigData(MatrixD & histData, int const i, int const j,
	             CopulaDef::CopInfoData * const  multiCopInf = nullptr,
	             DimT gridSize = 0);

	// at the moment, this is not supported
	// would probably have to add a private member with a matrix, since
	// matrix_row need a matrix as a parameter for the constructor
	//Cop2DData(VectorD const & marg1, VectorD const & marg2, int const nSamplPts,
	//          int const gridSize = 0);

	//void set_nmb_scens(DimT const nScens) override;

	/// rank-based cdf, using grid interpolation
	double cdfR(DimT const i, DimT const j) const override;

	virtual void set_nmb_scens(DimT const nScens) override;
};
*/


/*
class Cop2DData : public Cop2DComp {
private:
	/// poiner to the multivar info object
	CopulaDef::CopInfoData * p2multivarTg;
	DimT marg1idx;   ///< index of the first margin in the multivar info object
	DimT marg2idx;   ///< index of the first margin in the multivar info object

	MatrixF gridCdf; ///< the computed cdf values (values from [0,1])

	// have to have it, as it is pure virtual in the base - only a dummy
	double calc_cdf(double const u, double const v) const {
		throw std::logic_error("calc_cdf() should never be used in this class");
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

	void set_nmb_scens(DimT const nScens) override;
};
*/


/// Copula given by a 2D sample
/**
	\warning
	For the moment, this works only if all copulas are given by historical
	data, since it simply points to the multivariate \c CopInfoData object.
**//*
class Cop2DBigData : public Cop2DApproxGrid {
private:
	/// poiner to the multivar info object
	CopulaDef::CopInfoData * p2multivarTg;
	DimT marg1idx;   ///< index of the first margin in the multivar info object
	DimT marg2idx;   ///< index of the first margin in the multivar info object

	void calc_all_grid_cdfs() override; ///< overwriting the base method

protected:
	// - using matrix_row, as this does not allocate new data!
	// - matrix_row does not have a default (empty) constructor, so we cannot
	//   have an array - how would we initialize it?
	ublas::matrix_row<MatrixD> margin1; ///< values of the first margin
	ublas::matrix_row<MatrixD> margin2; ///< values of the first margin
	DimT nPts;     ///< number of the sample/data points in each margin


public:
	/// constructor with a matrix and row numbers
	Cop2DBigData(MatrixD & histData, int const i, int const j,
	             CopulaDef::CopInfoData * const  multiCopInf = nullptr,
	             DimT gridSize = 0);

	// at the moment, this is not supported
	// would probably have to add a private member with a matrix, since
	// matrix_row need a matrix as a parameter for the constructor
	//Cop2DData(VectorD const & marg1, VectorD const & marg2, int const nSamplPts,
	//          int const gridSize = 0);

	//void set_nmb_scens(DimT const nScens) override;
};
*/


/// Normal copula
/**
	\note I have not found any library providing the copula cdf, so I use the
	      bivariate distribution and then inverse of the marginal cdfs.
	      For the former, the only implementation with suitable license is from
	      QuantLib; for the latter, we use QuantLib, but boost have the function
	      as well...
**/
class Cop2DNormal : public Cop2DComp {
private:
	double correl;  ///< correlation

	/// \name QuantLib objects that compute normal CDFs etc.
	///@{
		/// inverse cdf (by default standard normal distrib.)
		QuantLib::InverseCumulativeNormal N01InvCdf;
		/// bivariate cdf (standardized, needs correlation)
		/**
			\note Must use a pointer, since this needs correlation on construction
			      and not all class constructors have this as a parameter.
		**/
		std::unique_ptr<QuantLib::BivariateCumulativeNormalDistribution> p2N01Cdf2D;
	///@}

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DNormal(double const rho);

	Cop2DNormal(std::istream & paramStr);
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
class Cop2DStudent : public Cop2DComp {
private:
	double correl;  ///< correlation
	unsigned dof;   ///< degree of freedom

	/// \name QuantLib objects that compute CDFs etc.
	/**
		\note Must use pointers, since the object needs dof on construction
		      and not all class constructors have this as a parameter.
	**/
	///@{
		///class providing bivariate student cdf
		std::unique_ptr<QuantLib::BivariateCumulativeStudentDistribution> p2tCdf2D;

		/// QuantLib object for the inverse cdf
		std::unique_ptr<QuantLib::InverseCumulativeStudent> p2tInvCdf;
	///@}

	double calc_cdf(double const u, double const v) const;

public:
	Cop2DStudent(unsigned degF, double const rho);

	Cop2DStudent(std::istream & paramStr);
};

} // namespace


#endif
