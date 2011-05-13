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
DimT Cop2DInfo::u_to_grid(double const u) const
{
	assert (gridN > 0 && "grid must be set up prior to calling u_to_grid()");
	assert (!customGridPts && "u_to_grid() does not work with custom grid pts");
	if (customGridPts) {
		// custom grid points -> would have to do a line search .. not done yet
		throw std::logic_error("call to u_to_grid() with custom grid points");
	}
	DimT i = (DimT) (u * gridN);
	assert (isEq(u, gridPts(i)) && "u must be a grid point");
	return i;
}

double Cop2DInfo::cdf(double const u, double const v) const
{
	if (useGrid) {
		return gridCdf(u_to_grid(u), u_to_grid(v));
	} else {
		assert (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && "bounds check");
		// handle boundaries manually - can be problematic in some cases!
		if (isEq(u, 0.0) || isEq(v, 0.0)) { return 0.0; }
		if (isEq(u, 1.0)) { return v; }
		if (isEq(v, 1.0)) { return u; }
		return calc_cdf(u,v);
	}
}


void Cop2DInfo::calc_all_grid_cdfs()
{
	assert (useGrid && gridN > 0 && gridPts.size() == gridN && "sanity checks");

	gridCdf.resize(gridN, gridN);
	for (DimT i = 0; i < gridN; ++i) {
		for (DimT j = 0; j < gridN; ++j) {
			gridCdf(i, j) = calc_cdf(gridPts(i), gridPts(j));
		}
	}
}


void Cop2DInfo::init_cdf_grid(DimT const N, double const posInInt)
{
	if (gridN == N && gridPts.size() == N && isEq(gridPts[0], posInInt / N)) {
		cerr << "Info: init_cdf_grid() called again with the same params!"
		     << endl;
		return; // the grid already exists, with the same specs.
	}
	// \todo bound checking

	useGrid = true;
	gridN = N;

	gridPts.resize(N);
	for (DimT i = 0; i < N; ++i) {
		gridPts[i] = (i + posInInt) / N;
		assert (u_to_grid(gridPts[i]) == i && "sanity check");
	}

	calc_all_grid_cdfs();
}


void Cop2DInfo::init_cdf_grid(UVector const & gridPos)
{
	if (gridN > 0)
		cerr << "Warning: calling init_cdf_grid() with existing grid!" << endl;
	// \todo bound checking

	gridN = gridPos.size();
	gridPts = gridPos;
	useGrid = true;
	customGridPts = true;

	gridCdf.resize(gridN, gridN);
	for (DimT i = 0; i < gridN; ++i) {
		for (DimT j = 0; j < gridN; ++j) {
			gridCdf(i, j) = calc_cdf(gridPts(i), gridPts(j));
		}
	}

	calc_all_grid_cdfs();
}


/*
void Cop2DInfo::attach_multivar_info(CopulaDef::CopInfoBy2D * const p2tg,
                                     DimT const i, DimT const j)
{
	p2multivarTg = p2tg;
	if (i >= 0) marg1idx = i;
	if (j >= 0) marg2idx = j;
}
*/


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

double Cop2DClayton::calc_cdf(const double u, double const v) const
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

double Cop2DMarshallOlkin::calc_cdf(const double u, const double v) const
{
	if (pow(u, alpha) >= pow(v, beta)) {
		return pow(u, 1.0 - alpha) * v;
	} else {
		return u * pow(v, 1.0 - beta);
	}
}


/*
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
*/


// -----------------------------------------------------------------------
// Copula given by a 2D sample
Cop2DData::Cop2DData(UMatrix & histData, int const i, int const j,
                     CopulaDef::CopInfoData * const p2CopInf)
: Cop2DInfo(),
  p2multivarTg(p2CopInf), marg1idx(i), marg2idx(j),
  margin1(histData, i), margin2(histData, j), nPts()
{
	/*
	if (gridN > 0) {
		init_grid();
	}
	*/
}


// has to overwrite the base method, since the values are computed recursively
void Cop2DData::calc_all_grid_cdfs()
{
	/// sample cdf evaluated on the grid; indexed -1 .. N-1
	/** implemented using boost::multi_array, to get indices starting from -1 **/
	IMatrix gridRCdf;
	typedef boost::multi_array_types::extent_range ExtRange;

	DimT i, j, s;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");
	// add -1 so we can use the recursive formula!
	gridRCdf.resize(boost::extents[ExtRange(-1,gridN)][ExtRange(-1,gridN)]);
	for (i = -1; i < gridN; i++) {
		for (j = -1; j < gridN; j++) {
			gridRCdf[i][j] = 0; // initialization .. should NOT be necessary!
		}
	}

	// get the ranks of the data
	if (p2multivarTg && marg1idx >= 0 && marg2idx >= 0) {
		// we have a pointer to the multivar info - ask for margin ranks there!
		ublas::matrix_row<UIMatrix> margRanks1
			(CopulaDef::cop_info_data_ranks(*p2multivarTg), marg1idx);
		ublas::matrix_row<UIMatrix> margRanks2
			(CopulaDef::cop_info_data_ranks(*p2multivarTg), marg2idx);

		// compute the grid-pdf, store it to gridRCdf
		for (s = 0; s < nPts; s++) {
			i = u012Rank( rank2U01(margRanks1(s), nPts), gridN );
			j = u012Rank( rank2U01(margRanks2(s), nPts), gridN );
			gridRCdf[i][j] ++;
		}
	} else {
		// we do not have access to ranks -> have to compute them first
		cerr << "rank computation not yet finished!" << endl; exit(1);
		/*
		UIVector margRanks1(nPts), margRanks2(nPts);
		get_ranks(margin1, margRanks1);
		get_ranks(margin2, margRanks2);

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
		*/
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

	// finally, put values to the main matrix
	gridCdf.resize(gridN, gridN);
	for (i = 0; i < gridN; ++i) {
		for (j = 0; j < gridN; ++j) {
			gridCdf(i, j) = (double) gridRCdf[i][j] / (double) gridN;
		}
	}
}

/*
double Cop2DData::cdf(const double u, const double v) const
{
	assert (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && "bounds check");
	if (isEq(u, 0.0) || isEq(v, 0.0)) { return 0.0; }
	if (isEq(u, 1.0)) { return v; }
	if (isEq(v, 1.0)) { return u; }
	assert (gridN > 0 && "Have to set grid size before calling cdf()!");
	int i = u012Rank(u, gridN);
	int j = u012Rank(v, gridN);
	assert ( i >= 0 && i < gridN && j >= 0 && j < gridN && "sanity check" );
	return (double) gridRCdf[i][j] / (double) nPts;
}
*/

// explicit instantiation - tell the compiler which instances to compile
//template class Cop2DData<Vector>;
//template class Cop2DData< std::vector<double> >; // the same as above


// -----------------------------------------------------------------------
// Normal copula
Cop2DNormal::Cop2DNormal(double const rho)
: Cop2DInfo(), correl(rho), N01Cdf2D(rho)
{
	if (rho < -1 || rho > 1) {
		throw std::out_of_range("correlation in normal copula out of range");
	}
}

double Cop2DNormal::calc_cdf(double const u, double const v) const
{
	double x = N01InvCdf(u);
	double y = N01InvCdf(v);
	return N01Cdf2D(x, y);
}

/*
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
*/

/*
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
*/

// explicit instantiation - tell the compiler which instances to compile
//template class Cop2DNormal<Vector>;
//template class Cop2DNormal< std::vector<double> >; // the same as above
