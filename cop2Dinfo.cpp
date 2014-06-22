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
/*
void Copula2D::make_2d_cop_name_map(Cop2DNameMapT_Old & cMap) {
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
*/


// -----------------------------------------------------------------------
// Cop2DInfo base class
// creates the map - new derived classes have to be added here!
void Cop2DInfo::init_name_map() {
	add_to_map<Cop2DIndep>("indep");
//	add_to_map<Cop2DData>("data");
//	add_to_map<Cop2DData>("sample");
	add_to_map<Cop2DNormal>("normal");
	add_to_map<Cop2DStudent>("student");
	add_to_map<Cop2DStudent>("t");
	add_to_map<Cop2DGumbel>("Gumbel");
	add_to_map<Cop2DClayton>("Clayton");
	add_to_map<Cop2DFrank>("Frank");
	add_to_map<Cop2DNelsen2>("Nelsen-2");
	add_to_map<Cop2DNelsen18>("Nelsen-18");
	add_to_map<Cop2DMarshallOlkin>("Marshall-Olkin");
	add_to_map<Cop2DMarshallOlkin>("M-O");
}

// template for adding new entries to the-name-map
template <class T>
void Cop2DInfo::add_to_map(string const & idStr)
{
	// throw an error if the element with the same id-string already exists
	if (nameMap.count(idStr) > 0) {
		throw std::logic_error("Tried to re-insert '" + idStr
		                       + "' into the name map.");
	}
	// convert to lower case
	string idStrLC = idStr;
	for (auto & ltr: idStrLC)
		ltr = std::tolower(ltr);
	// add to the map (using the new C++11 method)
	nameMap.emplace(idStrLC,
	                [](std::istream & paramStr) {
	                	return static_cast<Cop2DInfo*>(new T (paramStr));
	                });
}

// creates a new object based on its name
/// \todo checks if shapeName is in the map
Cop2DInfo * Cop2DInfo::make_new(string const & copName,
                                std::istream & paramStr)
{
	// create the map if it does not exist
	// this way, we do not need to create the map manually
	if (nameMap.size() == 0)
		init_name_map();
	// convert to lower case
	string copNameLC = copName;
	for (auto & ltr: copNameLC)
		ltr = std::tolower(ltr);
	if (nameMap.count(copNameLC) == 0) {
		string errMsg = "Unknown copula type '" + copName + "'. "
		                + "Supported copula ids are:";
		for (auto nm : nameMap)
			errMsg += " " + nm.first;
		throw std::runtime_error (errMsg + ".");
	}
	return Cop2DInfo::nameMap.at(copNameLC)(paramStr);
}


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

// initialize the map - required!
Cop2DNameMapT Cop2DInfo::nameMap;


// -----------------------------------------------------------------------
// Clayton copula

Cop2DClayton::Cop2DClayton(double const theta)
: Cop2DInfo(), th(theta)
{
	if (th < -1) {
		throw std::out_of_range("Clayton copula: parameter theta must be >= -1");
	}
}

Cop2DClayton::Cop2DClayton(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> th;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Clayton copula");
	}
	if (th < -1) {
		throw std::out_of_range("Clayton copula: parameter theta must be >= -1");
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
		throw std::out_of_range("Gumbel copula: parameter theta must be >= 1");
	}
}

Cop2DGumbel::Cop2DGumbel(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> th;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Gumbel copula");
	}
	if (th < 1) {
		throw std::out_of_range("Gumbel copula: parameter theta must be >= 1");
	}
	iTh = 1 / th;
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
		throw std::out_of_range("Frank copula: parameter theta must be >= -1");
	}
}

Cop2DFrank::Cop2DFrank(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> th;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Frank copula");
	}
	if (th < -1) {
		throw std::out_of_range("Frank copula: parameter theta must be >= -1");
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
		throw std::out_of_range("Nelsen-2 copula: parameter theta must be >= 1");
	}
}

Cop2DNelsen2::Cop2DNelsen2(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> th;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Nelsen-2 copula");
	}
	if (th < 1) {
		throw std::out_of_range("Nelsen-2 copula: parameter theta must be >= 1");
	}
}


// -----------------------------------------------------------------------
// copula 4.2.18 from Nelsen, pp. 118

Cop2DNelsen18::Cop2DNelsen18(double const theta)
: Cop2DInfo(), th(theta)
{
	if (th < 2) {
		throw std::out_of_range("Nelsen-18 copula: parameter theta must be >= 2");
	}
}

Cop2DNelsen18::Cop2DNelsen18(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> th;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Nelsen-18 copula");
	}
	if (th < 2) {
		throw std::out_of_range("Nelsen-18 copula: parameter theta must be >= 2");
	}
}


// -----------------------------------------------------------------------
// Marshall-Olkin copula from Nelsen, pp. 53
Cop2DMarshallOlkin::Cop2DMarshallOlkin(double const alpha_, double const beta_)
: Cop2DInfo(), alpha(alpha_), beta(beta_)
{
	if (alpha < 0 || alpha > 1) {
		throw std::out_of_range
			("Marshall-Olkin copula: param. alpha is out of range");
	}
	if (beta < 0 || beta > 1) {
		throw std::out_of_range
			("Marshall-Olkin copula: param. beta is out of range");
	}
}

Cop2DMarshallOlkin::Cop2DMarshallOlkin(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> alpha >> beta;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Marshall-Olkin copula");
	}
	if (alpha < 0 || alpha > 1) {
		throw std::out_of_range
			("Marshall-Olkin copula: param. alpha is out of range");
	}
	if (beta < 0 || beta > 1) {
		throw std::out_of_range
			("Marshall-Olkin copula: param. beta is out of range");
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
: Cop2DInfo(), correl(rho),
  p2N01Cdf2D(new QuantLib::BivariateCumulativeNormalDistribution(rho))
{
	if (rho < -1 || rho > 1) {
		throw std::out_of_range("correlation in normal copula out of range");
	}
}

Cop2DNormal::Cop2DNormal(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> correl;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the Gumbel copula");
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
  p2tCdf2D(new BivariateCumulativeStudentDistribution(degF, rho)),
  p2tInvCdf(new QuantLib::InverseCumulativeStudent(degF))
{
	if (dof == 0) {
		throw std::out_of_range("d.o.f. in student copula must be >= 1");
	}
	if (correl < -1 || correl > 1) {
		throw std::out_of_range("correlation in student copula out of range");
	}
}

Cop2DStudent::Cop2DStudent(istream & paramStr)
: Cop2DInfo()
{
	paramStr >> dof >> correl;
	if (! paramStr.good()) {
		throw std::runtime_error
			("error while reading parameters for the student copula");
	}
	if (dof == 0) {
		throw std::out_of_range("d.o.f. in student copula must be >= 1");
	}
	if (correl < -1 || correl > 1) {
		throw std::out_of_range("correlation in student copula out of range");
	}
	p2tCdf2D.reset(new BivariateCumulativeStudentDistribution(dof, correl));
	p2tInvCdf.reset(new QuantLib::InverseCumulativeStudent(dof));
}


double Cop2DStudent::calc_cdf(double const u, double const v) const
{
	double x = (*p2tInvCdf)(u);
	double y = (*p2tInvCdf)(v);
	return (*p2tCdf2D)(x, y);
}
