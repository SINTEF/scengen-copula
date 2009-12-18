#ifndef COP_2D_INFO_HPP
#define COP_2D_INFO_HPP

#include "common.hpp"

namespace Copula2D{

/// specifications of a bivariate copula, used as targets
class Cop2DInfo {
private:

protected:

public:
	Cop2DInfo() {}

	virtual ~Cop2DInfo() {}

	virtual double cdf(double const u, double const v) const = 0;
};


/// class for the transposed versions of all copulas
class Cop2DInfTr : public Cop2DInfo {
private:
	// note: this is the same as 'const Cop2DInfo *p2copInfo;'
	Cop2DInfo const *p2copInfo; ///< non-mutable pointer to the transposed obj.

public:
	Cop2DInfTr(Cop2DInfo const *const p2Transp) : p2copInfo(p2Transp) {}

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


/// copula 4.2.18 from Nelsen, pp. 118
class Cop2DNelsen18 : public Cop2DInfo {
private:
	double th; ///< parameter theta of the copula; th in [-1, inf)

public:
	Cop2DNelsen18(const double theta);

	virtual double cdf(const double u, const double v) const;
};


} // namespace

#endif
