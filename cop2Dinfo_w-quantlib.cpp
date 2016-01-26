/// implementation of 2D copulas that depend on QuantLib library
/**
	This makes it easier to create versions of the code that does neither
	include or depend on the QuantLib library

	For the moment, this file includes normal and student copulas
**/

#include "cop2Dinfo_w-quantlib.hpp"

#include <iostream>

//using namespace std;
using namespace Copula2D;


// -----------------------------------------------------------------------
// Normal copula

Cop2DNormal::Cop2DNormal(double const rho)
: Cop2DComp(), correl(rho),
  p2N01Cdf2D(new QuantLib::BivariateCumulativeNormalDistribution(rho))
{
	if (rho < -1 || rho > 1) {
		throw std::out_of_range("correlation in normal copula out of range");
	}
}

Cop2DNormal::Cop2DNormal(std::istream & paramStr)
: Cop2DComp()
{
	paramStr >> correl;
	if (paramStr.fail()) {
		throw std::runtime_error
			("error while reading parameters for the normal copula");
	}
	if (correl < -1 || correl > 1) {
		throw std::out_of_range("correlation in normal copula out of range");
	}
	p2N01Cdf2D.reset(new QuantLib::BivariateCumulativeNormalDistribution(correl));
}


double Cop2DNormal::calc_cdf(double const u, double const v) const
{
	double x = N01InvCdf(u);
	double y = N01InvCdf(v);
	return (*p2N01Cdf2D)(x, y);
}


// -----------------------------------------------------------------------
// Student-t copula

Cop2DStudent::Cop2DStudent(unsigned degF, double const rho)
: correl(rho), dof(degF),
  p2tCdf2D(new QuantLib::BivariateCumulativeStudentDistribution(degF, rho)),
  p2tInvCdf(new QuantLib::InverseCumulativeStudent(degF))
{
	if (dof == 0) {
		throw std::out_of_range("d.o.f. in student copula must be >= 1");
	}
	if (correl < -1 || correl > 1) {
		throw std::out_of_range("correlation in student copula out of range");
	}
}

Cop2DStudent::Cop2DStudent(std::istream & paramStr)
: Cop2DComp()
{
	paramStr >> dof >> correl;
	if (paramStr.fail()) {
		throw std::runtime_error
			("error while reading parameters for the student copula");
	}
	if (dof == 0) {
		throw std::out_of_range("d.o.f. in student copula must be >= 1");
	}
	if (correl < -1 || correl > 1) {
		throw std::out_of_range("correlation in student copula out of range");
	}
	p2tCdf2D.reset(new QuantLib::BivariateCumulativeStudentDistribution(dof, correl));
	p2tInvCdf.reset(new QuantLib::InverseCumulativeStudent(dof));
}


double Cop2DStudent::calc_cdf(double const u, double const v) const
{
	double x = (*p2tInvCdf)(u);
	double y = (*p2tInvCdf)(v);
	return (*p2tCdf2D)(x, y);
}
