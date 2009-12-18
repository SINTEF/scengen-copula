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
	return max(1 + th / log(exp(th / (u - 1)) + exp(th / (v - 1))), 0.0);
}
