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


// ---------------------------------------------------------------------------
// generic routines

// NB: Cop2DTypeID = {indep, Clayton, Gumbel, Frank, Nelsen2, Nelsen18,
//                    MarshallOlkin, data, normal, student};
void Copula2D::make_2d_cop_name_map(Cop2DNameMapT & cMap) {
	// note: when listing the map, it goes from last to first
	cMap["indep"]   = Cop2DTypeID::indep;
	cMap["Clayton"] = Cop2DTypeID::Clayton;
	cMap["clayton"] = Cop2DTypeID::Clayton;
	cMap["Gumbel"]  = Cop2DTypeID::Gumbel;
	cMap["gumbel"]  = Cop2DTypeID::Gumbel;
	cMap["Frank"]   = Cop2DTypeID::Frank;
	cMap["frank"]   = Cop2DTypeID::Frank;
	cMap["Nelsen-2"]       = Cop2DTypeID::Nelsen2;
	cMap["nelsen-2"]       = Cop2DTypeID::Nelsen2;
	cMap["Nelsen-18"]      = Cop2DTypeID::Nelsen18;
	cMap["nelsen-18"]      = Cop2DTypeID::Nelsen18;
	cMap["Marshall-Olkin"] = Cop2DTypeID::MarshallOlkin;
	cMap["marshall-olkin"] = Cop2DTypeID::MarshallOlkin;
	cMap["M-O"]     = Cop2DTypeID::MarshallOlkin;
	cMap["m-o"]     = Cop2DTypeID::MarshallOlkin;
	cMap["data"]    = Cop2DTypeID::data;
	cMap["sample"]  = Cop2DTypeID::data;
	cMap["normal"]  = Cop2DTypeID::normal;
	cMap["student"] = Cop2DTypeID::student;
	cMap["Student"] = Cop2DTypeID::student;
	cMap["t"]       = Cop2DTypeID::student;
}


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
	DimT i = (DimT) ((u - gridPts(0)) * gridN + DblEps);
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

double Cop2DInfo::cdfR(DimT const i, DimT const j) const
{
	if (useGrid) {
		return gridCdf(i, j);
	} else {
		// handle boundaries manually - can be problematic in some cases!
		if (i * j == 0) { return 0.0; }
		if (i == gridN - 1) { return gridPts(j); }
		if (j == gridN - 1) { return gridPts(i); }
		return calc_cdf(gridPts(i), gridPts(j));
	}
}


void Cop2DInfo::calc_all_grid_cdfs()
{
	assert (useGrid && "grid should be computed only if needed");
	assert (gridN > 0 && gridPts.size() == gridN && "sanity checks");

	useGrid = false; // needed to avoid cdf() trying to use the grid values
	gridCdf.resize(gridN, gridN);
	for (DimT i = 0; i < gridN; ++i) {
		for (DimT j = 0; j < gridN; ++j) {
			gridCdf(i, j) = cdf(gridPts(i), gridPts(j));
		}
	}
	useGrid = true; // reset back
}


//! gridPts is used also when we do not use grid-based cdf!
//! \todo rename useGrid to usedGridCdf? !!
void Cop2DInfo::init_cdf_grid(DimT const N, double const posInInt)
{
	if (gridN == N && gridPts.size() == N && isEq(gridPts[0], posInInt / N)) {
		// this happens with independent copula, which re-uses the 2D object
		TRACE (TrInfo, "init_cdf_grid() called again with the same params");
		return; // the grid already exists, with the same specs.
	}
	// \todo bound checking

	gridN = N;

	gridPts.resize(N);
	for (DimT i = 0; i < N; ++i) {
		gridPts(i) = (i + posInInt) / N;
		assert (u_to_grid(gridPts(i)) == i && "sanity check");
	}
	TRACE (TrDetail3, "posInInt = " << posInInt << "; gridPts = " << gridPts);

	if (useGrid) {
		calc_all_grid_cdfs();
	}
}


void Cop2DInfo::init_cdf_grid(VectorD const & gridPos)
{
	if (! useGrid) {
		return;
	}
	if (gridN > 0)
		cerr << "Warning: calling init_cdf_grid() with existing grid!" << endl;
	// \todo bound checking

	gridN = gridPos.size();
	gridPts = gridPos;
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
void Cop2DInfo::clear_cdf_grid()
{
	gridCdf.resize(0, 0, false); // clear() does not free memory?!
	useGrid = false;
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
// Gumbel copula

Cop2DGumbel::Cop2DGumbel(double const theta)
: Cop2DInfo(), th(theta), iTh(1 / theta)
{
	if (th < 1) {
		cerr << "ERROR: Parameter theta out of range!";
		exit(1);
	}
}

double Cop2DGumbel::calc_cdf(const double u, double const v) const
{
	return exp(-pow(pow(-log(u), th) + pow(-log(v), th), iTh));
}


// -----------------------------------------------------------------------
// Frank copula

Cop2DFrank::Cop2DFrank(double const theta)
: Cop2DInfo(), th(theta), C(exp(-theta) - 1)
{
	if (th < -1) {
		cerr << "ERROR: Parameter theta out of range!";
		exit(1);
	}
}

double Cop2DFrank::calc_cdf(const double u, double const v) const
{
	if (fabs(th) < DblEps) {
		// Frank cdf is undefined for theta = 0
		// but the limit of the cdf for theta -> 0 is the indep. copula
		return u * v;
	} else {
		return -1 * log(1 + (exp(-th * u) - 1) * (exp(-th * v) - 1) / C) / th;
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


// -----------------------------------------------------------------------
// Copula given by a 2D sample
Cop2DData::Cop2DData(MatrixD & histData, int const i, int const j,
                     CopulaDef::CopInfoData * const p2CopInf)
: Cop2DInfo(),
  p2multivarTg(p2CopInf), marg1idx(i), marg2idx(j),
  margin1(histData, i), margin2(histData, j), nPts(histData.size2())
{}


// has to overwrite the base method, since the values are computed recursively
void Cop2DData::calc_all_grid_cdfs()
{
	int i, j;
	DimT s;
	assert (gridN > 0 && "Have to set grid size before calling init_grid()!");

	/// sample cdf evaluated on the grid; indexed -1 .. N-1
	/** implemented using boost::multi_array, to get indices starting from -1 **/
	Array2D<int> gridRCdf(gridN, gridN, -1);
	//boost::multi_array<int, 2> gridRCdf;
	//typedef boost::multi_array_types::extent_range ExtRange;

	// add -1 so we can use the recursive formula!
	// gridRCdf.resize(boost::extents[ExtRange(-1,gridN)][ExtRange(-1,gridN)]);
	for (i = -1; i < (int) gridN; i++) {
		for (j = -1; j < (int) gridN; j++) {
			gridRCdf(i,j) = 0; // initialization .. should NOT be necessary!
		}
	}

	// get the ranks of the data
	if (p2multivarTg && marg1idx >= 0 && marg2idx >= 0) {
		// we have a pointer to the multivar info - ask for margin ranks there!
		ublas::matrix_row<MatrixI> margRanks1
			(CopulaDef::cop_info_data_ranks(*p2multivarTg), marg1idx);
		ublas::matrix_row<MatrixI> margRanks2
			(CopulaDef::cop_info_data_ranks(*p2multivarTg), marg2idx);

		// compute the grid-pdf, store it to gridRCdf
		for (s = 0; s < nPts; s++) {
			i = u012Rank( rank2U01(margRanks1(s), nPts), gridN );
			j = u012Rank( rank2U01(margRanks2(s), nPts), gridN );
			gridRCdf(i,j) ++;
		}
	} else {
		// we do not have access to ranks -> have to compute them first
		cerr << "rank computation not yet finished!" << endl; exit(1);
		/*
		VectorI margRanks1(nPts), margRanks2(nPts);
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
	for (i = 0; i < (int) gridN; i++) {
		int colCount = 0;
		for (j = 0; j < (int) gridN; j++) {
			// at this point, gridRCdf(i,j) includes the rank pdf from above
			colCount += gridRCdf(i,j);
			gridRCdf(i,j) = gridRCdf(i-1,j) + colCount;
			// at this point, gridRCdf(i,j) includes the rank cdf!
		}
	}

	// finally, put values to the main matrix
	gridCdf.resize(gridN, gridN);
	for (i = 0; i < (int) gridN; ++i) {
		for (j = 0; j < (int) gridN; ++j) {
			gridCdf(i, j) = (double) gridRCdf(i,j) / (double) nPts; // was: gridN;
		}
	}
}


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


// -----------------------------------------------------------------------
// Student-t copula
Cop2DStudent::Cop2DStudent(unsigned degF, double const rho)
: correl(rho), dof(degF), tCdf2D(degF, rho), tInvCdf(degF)
{}


double Cop2DStudent::calc_cdf(double const u, double const v) const
{
	double x = tInvCdf(u);
	double y = tInvCdf(v);
	return tCdf2D(x, y);
}
