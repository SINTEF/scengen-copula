#ifndef COPULA_INFO_HPP
#define COPULA_INFO_HPP

#include "common.hpp"
#include "cassert"

namespace CopulaDef{

/// specifications of a bivariate copula, used as targets
class CopulaInfo {
private:

protected:

public:
	CopulaInfo() {}

	virtual ~CopulaInfo() {}

	/// cdf at vector \a u = (u_1, .. , u_n)
	virtual double cdf(TVectorD const u) const = 0;
};


/// terget copula described by a historical data (sample)
class CopInfoData : public CopulaInfo {
private:
	TMatrixI & ranks; ///< reference to the matrix of margin ranks
	TMatrixD u01Data; ///< the ranks transformed to U(0,1)

protected:
	int nVars; ///< number of random variables
	int nPts; ///< number of the sample/data points

public:
	CopInfoData(TMatrixI & ranksMat, bool const fillU01Data = true);

	~CopInfoData() {}

	double cdf(TVectorD const u) const;
};

} // namespace

#endif
