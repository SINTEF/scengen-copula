#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <algorithm>
//#include <deque>
//#include <vector>
//#include <cassert>

#include "cop2Dinfo.hpp"

using namespace std;
using namespace Copula2D;


// -------------------------------------------------------------
// Cop2DInfo base class


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


// -------------------------------------------------------------
// Copula given by a 2D sample
template <class T>
Cop2DData<T>::Cop2DData(T const * marg1, T const * marg2, int const nSamplPts,
												int const gridSize)
: Cop2DInfo(),
  gridN(gridSize), nPts(nSamplPts)
{
	margins[0] = marg1;
	margins[1] = marg2;
	if (gridN > 0) {
		init_grid();
	}
}

template <class T>
int Cop2DData<T>::init_grid()
{
	int i, j, s;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");
	// add -1 so we can use the recursive formula!
	gridRCdf.resize(boost::extents[Range(-1,gridN)][Range(-1,gridN)]);
	for (i = -1; i < gridN; i++) {
		for (j = -1; j < gridN; j++) {
			gridRCdf[i][j] = 0; // initialization .. should NOT be necessary!
		}
	}

	// get the vectors of ranks
	std::vector<int> margRanks[2];
	for (int marg = 0; marg < 2; marg++) {
		//std::vector<int> ranks;
		margRanks[marg].resize(nPts);
		get_ranks(*(margins[marg]), margRanks[marg]);
	}

	// compute the grid-pdf, store it to gridRCdf
	double tmpF = (double) gridN / (double) nPts;
	for (s = 0; s < nPts; s++) {
		i = floor((margRanks[0][s]) * tmpF);
		j = floor((margRanks[1][s]) * tmpF);
		gridRCdf[i][j] ++;
	}

	// compute the grid-cdf
	for (i = 0; i < gridN; i++) {
		int colCount = 0;
		for (j = 0; j < gridN; j++) {
			// at this point, gridRCdf[i][j] includes the rank pdf from above
			colCount += gridRCdf[i][j];
			gridRCdf[i][j] = gridRCdf[i-1][j] + colCount;
			// at this point, gridRCdf[i][j] includes the rank cdf!
		}
	}

	return 0;
}

template <class T>
double Cop2DData<T>::cdf(const double u, const double v) const
{
	int i = min((int) floor(u * gridN), gridN - 1);
	int j = min((int) floor(v * gridN), gridN - 1);
	//cerr << "TMP: Cop2DData::cdf(): i = " << i << "; j = " << j << endl;
	//cerr << "TMP: gridRCdf[i][j] = " << gridRCdf[i][j] << endl;
	return (double) gridRCdf[i][j] / (double) nPts;
}

// explicit instantiation - tell the compiler which instances to compile
template class Cop2DData<Vector>;
//template class Cop2DData< std::vector<double> >; // the same as above
