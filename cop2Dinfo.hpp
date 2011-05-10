#ifndef COP_2D_INFO_HPP
#define COP_2D_INFO_HPP

#include <cassert>
// QuantLib libraries used for the normal copula
#include <ql/math/distributions/bivariatenormaldistribution.hpp>

#include "common.hpp"


// forward declarations, needed by the following class
namespace CopulaDef {
	class CopInfoBy2D;
	class CopInfoData;
	UMatrix & cop_info_data_vals(CopInfoData & copInfo);
	UIMatrix & cop_info_data_ranks(CopInfoData & copInfo);
}

namespace Copula2D{

/// specifications of a bivariate copula, used as targets
class Cop2DInfo {
private:

protected:
	/// poiner to the multivar info object
	CopulaDef::CopInfoBy2D * p2multivarTg;
	DimT marg1idx; ///< index of the first margin in the multivar info object
	DimT marg2idx; ///< index of the first margin in the multivar info object

public:
	Cop2DInfo(CopulaDef::CopInfoBy2D * const p2CopInfo = NULL,
	          DimT const i = -1, DimT const j = -1)
	: p2multivarTg(p2CopInfo), marg1idx(i), marg2idx(j)
	{}

	virtual ~Cop2DInfo() {}

	/// cdf - calculated, so slower than the grid-based version
	virtual double cdf(double const u, double const v) const = 0;

	typedef boost::shared_ptr<Cop2DInfo> Ptr;

	void attach_multivar_info(CopulaDef::CopInfoBy2D * const p2tg,
	                          DimT const i = -1, DimT const j = -1)
	{
		p2multivarTg = p2tg;
		if (i >= 0) marg1idx = i;
		if (j >= 0) marg2idx = j;
	}
};


/// class for the transposed versions of all copulas
class Cop2DInfTr : public Cop2DInfo {
private:
	// note: this is the same as 'const Cop2DInfo *p2copInfo;'
	Cop2DInfo const *p2copInfo; ///< non-mutable pointer to the transposed obj.

public:
	/// constructor using a (raw) pointer to the transposed object
	Cop2DInfTr(Cop2DInfo const & p2Transp)
	: Cop2DInfo(p2Transp), p2copInfo(&p2Transp) {}

	/// constructor using a smart pointer to the transposed object
	// the .get() method returns a raw pointer from a smart one
	Cop2DInfTr(Cop2DInfo::Ptr const p2Transp) : p2copInfo(p2Transp.get()) {}

	double cdf(double const u, double const v) const {
		return p2copInfo->cdf(v,u);
	}
};


//-----------------------------------------------------------------------
//  DERIVED CLASSES

/// The independent copula
class Cop2DIndep : public Cop2DInfo {
public:
	Cop2DIndep() : Cop2DInfo() {}

	double cdf(const double u, const double v) const { return u * v; }
};


/// The Clayton copula
/// copula 4.2.1 from Nelsen, pp. 116
class Cop2DClayton : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf)

public:
	Cop2DClayton(const double theta);

	virtual double cdf(const double u, const double v) const;
};


/// copula 4.2.2 from Nelsen, pp. 116
class Cop2DNelsen2 : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [1, inf)

public:
	Cop2DNelsen2(const double theta);

	virtual double cdf(const double u, const double v) const;
};


/// copula 4.2.18 from Nelsen, pp. 118
class Cop2DNelsen18 : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf)

public:
	Cop2DNelsen18(const double theta);

	virtual double cdf(const double u, const double v) const;
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

public:
	Cop2DMarshallOlkin(double const alpha_, double const beta_);

	virtual double cdf(const double u, const double v) const;
};


/// Copula given by a 2D sample
/**
	T is the type of the margin vector (double * or std::vector<double>)
**/
template <class T>
class Cop2DDataOld : public Cop2DInfo {
private:
	int gridN;        ///< size of the grid on which we compute the pdf and cdf
	IMatrix gridRCdf; ///< sample cdf evaluated on the grid; indexed -1 .. N-1

protected:
	T const * margins[2];  ///< the two vectors of margins
	int nPts; ///< number of the sample/data points

	int init_grid();

public:
	Cop2DDataOld(T const * marg1, T const * marg2, int const nSamplPts,
	          int const gridSize = 0);

	int set_grid_size(int const N) { gridN = N; return init_grid(); }

	virtual double cdf(const double u, const double v) const;
};


/// Copula given by a 2D sample
class Cop2DData : public Cop2DInfo {
private:
	CopulaDef::CopInfoData * p2multivarTg;
	int gridN; ///< size of the grid on which we compute the pdf and cdf
	/// sample cdf evaluated on the grid; indexed -1 .. N-1
	/** implemented using boost::multi_array, to get indices starting from -1 **/
	IMatrix gridRCdf;

protected:
	// - using matrix_row, as this does not allocate new data!
	// - matrix_row does not have a default (empty) constructor, so we cannot
	//   have an array - how would we initialize it?
	ublas::matrix_row<UMatrix> margin1; ///< values of the first margin
	ublas::matrix_row<UMatrix> margin2; ///< values of the first margin
	int nPts;      ///< number of the sample/data points in each margin

	void init_grid();

public:
	/// constructor with a matrix and row numbers
	Cop2DData(UMatrix & histData, int const i, int const j,
	          int const gridSize = 0,
	          CopulaDef::CopInfoData * const  multiCopInf = NULL);

	// at the moment, this is not supported
	// would probably have to add a private member with a matrix, since
	// matrix_row need a matrix as a parameter for the constructor
	//Cop2DData(UVector const & marg1, UVector const & marg2, int const nSamplPts,
	//          int const gridSize = 0);

	void set_grid_size(int const N) { gridN = N; init_grid(); }

	virtual double cdf(const double u, const double v) const;
};


/// Normal copula
/**
	T is the type of the margin vector (double * or std::vector<double>)
**/
template <class T>
class Cop2DNormal : public Cop2DInfo {
private:
	double correl;  ///< correlation
	int gridN;      ///< size of the grid on which we compute the pdf and cdf
	Matrix gridCdf; ///< cdf evaluated on the grid; indexed -1 .. N-1

	/// \name QuantLib objects that compute normal CDFs etc.
	///@{
		/// inverse cdf (by default standard normal distrib.)
		QuantLib::InverseCumulativeNormal N01InvCdf;
		/// bivariate cdf (standardized, needs correlation)
		QuantLib::BivariateCumulativeNormalDistribution N01Cdf2D;
	///@}

protected:
	int init_grid();

public:
	Cop2DNormal(double const rho, int const gridSize = 0);

	int set_grid_size(int const N) { gridN = N; return init_grid(); }

	/// computes the cdf of a given pair of values
	/**
		If \c gridN > 0, it uses the values from the grid, otherwise it
		computes the cdf directly - which is probably slower.

		\warning The two variants give (obviously) the same results only on
		         the grid. This is important as the grid is in the borders of
		         the intervals of the support points, but the \a u and \a v
		         values are (by default) in the middle .. so we never get 0 or 1,
		         which we can get on the grid! <br>
		         All in all, it is probably better to use the grid...
	**/
	virtual double cdf(const double u, const double v) const;
};


} // namespace

#endif
