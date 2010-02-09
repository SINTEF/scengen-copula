#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>

#include "cop2Dinfo.hpp"

using namespace std;
using namespace Copula2D;

// -------------------------------------------------------------
// Clayton copula

Cop2DClayton::Cop2DClayton(double const theta)
: Cop2DInfo(), th(theta)
{
	if (th < -1) {
		cerr << "ERROR: Parameter theta out of range!";
		exit(1);
	}
}

double Cop2DClayton::cdf(const double u, double const v) const
{
	if (fabs(th) < DblEps) {
		// Clayton cdf is undefined for theta = 0
		// but the limit of the cdf for theta -> 0 is the indep. copula
		return u * v;
	} else {
		return pow(max(pow(u, -th) + pow(v, -th) - 1.0, 0.0), -1.0 / th);
	}
}


// -------------------------------------------------------------
// copula 4.2.18 from Nelsen, pp. 118

Cop2DNelsen18::Cop2DNelsen18(double const theta)
: Cop2DInfo(), th(theta)
{
	if (th < 2) {
		cerr << "ERROR: Parameter theta out of range!";
		exit(1);
	}
}

double Cop2DNelsen18::cdf(const double u, double const v) const
{
	// avoiding numerical issues
	if (u < DblEps || v < DblEps)
		return 0.0;
	if (u + DblEps > 1.0)
		return v;
	if (v + DblEps > 1.0)
		return u;

	return max(1 + th / log(exp(th / (u - 1)) + exp(th / (v - 1))), 0.0);
}


// -------------------------------------------------------------
// Marshall-Olkin copula from Nelsen, pp. 53
Cop2DMarshallOlkin::Cop2DMarshallOlkin(double const alpha_, double const beta_)
: Cop2DInfo(), alpha(alpha_), beta(beta_)
{
	if (alpha < 0 || alpha > 1) {
		cerr << "ERROR: Parameter alpha out of range!";
		exit(1);
	}
	if (beta < 0 || beta > 1) {
		cerr << "ERROR: Parameter beta out of range!";
		exit(1);
	}
}

double Cop2DMarshallOlkin::cdf(const double u, const double v) const
{
	if (pow(u, alpha) >= pow(v, beta)) {
		return pow(u, 1.0 - alpha) * v;
	} else {
		return u * pow(v, 1.0 - beta);
	}
}

