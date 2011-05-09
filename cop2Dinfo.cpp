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


// -----------------------------------------------------------------------
// Cop2DInfo base class


// -----------------------------------------------------------------------
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


// -----------------------------------------------------------------------
// copula 4.2.2 from Nelsen, pp. 116

Cop2DNelsen2::Cop2DNelsen2(double const theta)
: Cop2DInfo(), th(theta)
{
	if (th < 1) {
		cerr << "ERROR: Parameter theta out of range!";
		exit(1);
	}
}

double Cop2DNelsen2::cdf(const double u, double const v) const
{
	// avoiding numerical issues
	if (u < DblEps || v < DblEps)
		return 0.0;
	if (u + DblEps > 1.0)
		return v;
	if (v + DblEps > 1.0)
		return u;

	return max(1 - pow(pow(1 - u, th) + pow(1 - v, th), 1.0 / th), 0.0);
}


// -----------------------------------------------------------------------
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


// -----------------------------------------------------------------------
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


// -----------------------------------------------------------------------
// Copula given by a 2D sample
template <class T>
Cop2DDataOld<T>::Cop2DDataOld(T const * marg1, T const * marg2, int const nSamplPts,
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
int Cop2DDataOld<T>::init_grid()
{
	int i, j, s;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");
	// add -1 so we can use the recursive formula!
	gridRCdf.resize(boost::extents[ExtRange(-1,gridN)][ExtRange(-1,gridN)]);
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
		/// \todo Use the already computed ranks!
		get_ranks(*(margins[marg]), margRanks[marg]);
	}

	// compute the grid-pdf, store it to gridRCdf
	for (s = 0; s < nPts; s++) {
		i = u012Rank( rank2U01(margRanks[0][s], nPts), gridN );
		j = u012Rank( rank2U01(margRanks[1][s], nPts), gridN );
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
double Cop2DDataOld<T>::cdf(const double u, const double v) const
{
	assert (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && "bounds check");
	if (isEq(u, 0.0) || isEq(v, 0.0)) { return 0.0; }
	if (isEq(u, 1.0)) { return v; }
	if (isEq(v, 1.0)) { return u; }
	int i = u012Rank(u, gridN);
	int j = u012Rank(v, gridN);
	assert ( i >= 0 && i < gridN && j >= 0 && j < gridN && "sanity check" );
	return (double) gridRCdf[i][j] / (double) nPts;
}

// explicit instantiation - tell the compiler which instances to compile
template class Cop2DDataOld<Vector>;
//template class Cop2DDataOld< std::vector<double> >; // the same as above


// -----------------------------------------------------------------------
// Copula given by a 2D sample
Cop2DData::Cop2DData(UVector const & marg1, UVector const & marg2,
                     int const nSamplPts, int const gridSize)
: Cop2DInfo(),
  gridN(gridSize), margin1(marg1), margin2(marg2), nPts(nSamplPts)
{
	if (gridN > 0) {
		init_grid();
	}
}

int Cop2DData::init_grid()
{
	int i, j, s;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");
	// add -1 so we can use the recursive formula!
	gridRCdf.resize(boost::extents[ExtRange(-1,gridN)][ExtRange(-1,gridN)]);
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
		/// \todo Use the already computed ranks!
		get_ranks(*(margins[marg]), margRanks[marg]);
	}

	// compute the grid-pdf, store it to gridRCdf
	for (s = 0; s < nPts; s++) {
		i = u012Rank( rank2U01(margRanks[0][s], nPts), gridN );
		j = u012Rank( rank2U01(margRanks[1][s], nPts), gridN );
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
	assert (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && "bounds check");
	if (isEq(u, 0.0) || isEq(v, 0.0)) { return 0.0; }
	if (isEq(u, 1.0)) { return v; }
	if (isEq(v, 1.0)) { return u; }
	int i = u012Rank(u, gridN);
	int j = u012Rank(v, gridN);
	assert ( i >= 0 && i < gridN && j >= 0 && j < gridN && "sanity check" );
	return (double) gridRCdf[i][j] / (double) nPts;
}

// explicit instantiation - tell the compiler which instances to compile
template class Cop2DData<Vector>;
//template class Cop2DData< std::vector<double> >; // the same as above


// -----------------------------------------------------------------------
// Normal copula
template <class T>
Cop2DNormal<T>::Cop2DNormal(double const rho, int const gridSize)
: Cop2DInfo(), correl(rho), gridN(gridSize), N01Cdf2D(rho)
{
	if (gridN > 0) {
		init_grid();
	}
}

template <class T>
int Cop2DNormal<T>::init_grid()
{
	int i, j;
	double u, v, x, y;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");
	// add -1 so we can use the recursive formula!
	gridCdf.resize(boost::extents[ExtRange(-1,gridN)][ExtRange(-1,gridN)]);
	for (i = -1; i < gridN; i++) {
		for (j = -1; j < gridN; j++) {
			gridCdf[i][j] = 0; // initialization .. should NOT be necessary!
		}
	}
	for (i = -1; i < gridN; ++i) {
		gridCdf[i][-1] = 0.0;
		gridCdf[-1][i] = 0.0;
		gridCdf[i][gridN-1] = static_cast<double>(i + 1) / gridN;
		gridCdf[gridN-1][i] = static_cast<double>(i + 1) / gridN;
	}

	for (i = 0; i < gridN - 1; i++) {
		u = rank2U01(i, gridN);
		x = N01InvCdf(u);
		for (j = 0; j < gridN - 1; j++) {
			v = rank2U01(j, gridN);
			y = N01InvCdf(v);
			gridCdf[i][j] = N01Cdf2D(x, y);
		}
	}

	return 0;
}

template <class T>
double Cop2DNormal<T>::cdf(const double u, const double v) const
{
	double F;
	assert (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && "bounds check");
	if (isEq(u, 0.0) || isEq(v, 0.0)) { return 0.0; }
	if (isEq(u, 1.0)) { return v; }
	if (isEq(v, 1.0)) { return u; }
	if (gridN > 0) {
		// there is a grid -> use it
		int i = u012Rank(u, gridN);
		int j = u012Rank(v, gridN);
		assert ( i >= 0 && i < gridN && j >= 0 && j < gridN && "sanity check" );
		F = gridCdf[i][j];
	} else {
		// no grid -> compute it (probably slower)
		double x = N01InvCdf(u);
		double y = N01InvCdf(v);
		F = N01Cdf2D(x, y);
	}
	//cout << "TMP: F(" << u << ", " << v << ") = " << F << endl;
	return F;
}

// explicit instantiation - tell the compiler which instances to compile
template class Cop2DNormal<Vector>;
//template class Cop2DNormal< std::vector<double> >; // the same as above
