#include "common.hpp"
#include "copula-info.hpp"
#include "cop2Dinfo.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace CopulaDef;
using std::cout;
using std::cerr;
using std::endl;


// ---------------------------------------------------------------------------
// generic routines

void CopulaDef::make_cop_name_map(CopNameMapT & cMap) {
	// note: when listing the map, it goes from last to first
	cMap["sample"] = CopTypeID::sample;
	cMap["normal"] = CopTypeID::normal;
	cMap["indep"] = CopTypeID::indep;
	cMap["student"] = CopTypeID::student;
	cMap["t"] = CopTypeID::student;
	cMap["mixed"] = CopTypeID::mixed;
}


// ---------------------------------------------------------------------------
// class CopulaInfo
void CopulaInfo::get_correl_matrix_from_stream(
	std::istream & is, ublas::symmetric_matrix<double> & X)
{
	is >> X; // this should read the matrix, including dimensions

	nVars = X.size1();
	assert (X.size2() == nVars && "must be a square matrix");

	#ifndef NDEBUG
		for (DimT i = 0; i < nVars; ++i) {
			if (! isEq(X(i, i), 1.0)) {
				throw std::range_error("correl(i,i) must be 1");
			}
			for (DimT j = 0; j < i; ++j) {
				double rho = X(i, j);
				if (rho < -1 || rho > 1 || !isEq(rho, X(j, i))) {
					throw std::range_error("error in the correlation matrix");
				}
			}
		}
	#endif
}


// ---------------------------------------------------------------------------
// class CopInfoIndep

double CopInfoIndep::cdf(VectorD const u) const {
	assert (u.size() == nVars && "dimension check");

	double F = u[0];
	DimT i;
	for (i = 1; i < nVars; ++i)
		F *= u[i];

	return F;
}


void CopInfoIndep::setup_2d_targets()
{
	assert (nVars > 0 && "must have a known dimension by now");
	if (p2Info2D.size1() != nVars || p2Info2D.size2() != nVars) {
		if (p2Info2D.size1() * p2Info2D.size2() > 0) {
			cerr << "Warning: resizing the matrix of 2D-copula objects!" << endl;
		}
		p2Info2D.resize(nVars, nVars);
	}

	p2bivarIndep = boost::make_shared<Copula2D::Cop2DIndep>();

	DimT i, j;
	for (i = 0; i < nVars; ++i) {
		for (j = i + 1; j < nVars; ++j) {
			p2Info2D(i, j) = p2bivarIndep;
			n2Dcops++;
		}
	}
}


// constructor with dimension
CopInfoIndep::CopInfoIndep(DimT const N)
: CopInfoBy2D(N, true), p2bivarIndep(new Copula2D::Cop2DIndep())
{
	setup_2d_targets();
}

// constructor with file name of the target distribution
CopInfoIndep::CopInfoIndep(std::string const & tgFName)
: CopInfoBy2D(0, true)
{
	// read the input file
	std::ifstream tgFStr(tgFName.c_str());
	if (!tgFStr) {
		throw std::ios_base::failure("Could not open input file `"
		                             + tgFName + "'!");
	}
	tgFStr >> nVars; // the only value we need

	setup_2d_targets();
}


// ---------------------------------------------------------------------------
// class CopInfoBy2D

/*
// initialize cdf grids for all the target 2D copulas; regular intervals
void CopInfoBy2D::init_cdf_grids(DimT const N, bool const useTgPos,
                                 double const posInInt)
{
	DimT i, j;
	for (i = 0; i < nVars; ++i) {
		for (j = i+1; j < nVars; ++j) {
			if (p2Info2D(i, j)) {
				assert (p2Info2D(i, j).get() != nullptr && "testing the null-check");
				if (useTgPos)
					p2Info2D(i, j)->init_cdf_grid(N);
				else
					p2Info2D(i, j)->init_cdf_grid(N, posInInt);
			}
		}
	}
}

// initialize cdf grids for all the target 2D copulas; custom grid points
void CopInfoBy2D::init_cdf_grids(VectorD const & gridPos)
{
	DimT i, j;
	for (i = 0; i < nVars; ++i) {
		for (j = i+1; j < nVars; ++j) {
			if (p2Info2D(i, j)) {
				assert (p2Info2D(i, j).get() != nullptr && "testing the null-check");
				p2Info2D(i, j)->init_cdf_grid(gridPos);
			}
		}
	}
}
*/


// set the number of scenarios (passed to the bivariate copulas)
void CopInfoBy2D::set_nmb_scens(DimT const nScens)
{
	DimT i, j;
	for (i = 0; i < nVars; ++i) {
		for (j = i+1; j < nVars; ++j) {
			if (p2Info2D(i, j)) {
				p2Info2D(i, j)->set_nmb_scens(nScens);
			}
		}
	}
}


// get the number of initialized 2D target copulas
DimT CopInfoBy2D::get_nmb_2d_copulas()
{
	DimT i, j;
	DimT nmbC = 0;

	for (i = 0; i < nVars; ++i) {
		for (j = i+1; j < nVars; ++j) { // p2Info2D is an upper-triang. matrix
			if (p2Info2D(i, j)) {
				nmbC++;
			}
		}
	}
	return nmbC;
}


// read the list of copula pairs from a file
std::list<std::pair<DimT, DimT>> CopInfoBy2D::read_pair_list(std::string const & fName)
{
	std::ifstream inFStream(fName.c_str());
	if (!inFStream) {
		throw std::ios_base::failure("Could not open input file `" + fName + "'!");
	}
	std::list<std::pair<DimT, DimT>> list;
	std::string line;
	while (std::getline(inFStream, line)) {
		if (line[0] == '#') {
			TRACE(TrInfo3, "skipping comment line in the copula-pairs file");
			continue;
		}
		std::stringstream lineStream(line);
		DimT i, j;
		lineStream >> i;
		lineStream >> j;
		// NB: margins should be indexed from 1 in the file
		if (i < 1 || j < 1) {
			throw std::domain_error("margin numbers in the copula pairs list must be indexed from 1!");
		}
		i--;
		j--;
		if (i == j) {
			throw std::invalid_argument("cannot create a copula of a margin with itself!");
		}
		if (i < j) {
			list.emplace_back(i, j);
		}
		else {
			list.emplace_back(j, i);
		}
	}
	
	return list;
}


// ---------------------------------------------------------------------------
// class CopInfoGen2D

CopInfoGen2D::CopInfoGen2D(std::string const & tgFName)
: CopInfoBy2D(0, false)
{
	// read the input file
	std::ifstream tgFStr(tgFName.c_str());
	if (!tgFStr) {
		throw std::ios_base::failure("Could not open input file `"
		                             + tgFName + "'!");
	}
	tgFStr >> nVars; // the first parameter is the dimension
	p2Info2D.resize(nVars, nVars);

	//Copula2D::Cop2DNameMapT_Old copNameMap;  ///< convert copula name to type
	//make_2d_cop_name_map(copNameMap);
	Copula2D::Cop2DInfo::Ptr p2tgCop;

	while (! tgFStr.eof()) {
		DimT i, j;
		tgFStr >> i >> j;
		if (! tgFStr.good()) {
			if (n2Dcops == 0)
				throw std::runtime_error("data file for mixed copulas should"
				                         "start with 3 numbers: dimension and "
				                         "the margins of the first 2D-copula");
			else
				break; // stop reading
		}
		if (i < 1 || i > nVars)
			throw std::runtime_error("invalid 1st margin for 2D-copula no. "
			                         + std::to_string(n2Dcops + 1));
		if (j < 1 || j > nVars)
			throw std::runtime_error("invalid 2nd margin for 2D-copula no. "
			                         + std::to_string(n2Dcops + 1));
		--i; --j; // assume the inputs count from 1, not from zero
		std::string copType;
		tgFStr >> copType;
		if (! tgFStr.good())
			throw std::runtime_error("wrong input format for 2D-copula "
			                         + std::to_string(n2Dcops + 1)
			                         + ": expected copula type (string)");

		// NEW - NOT YET FINISHED!
		std::string paramsAsString;
		std::getline(tgFStr, paramsAsString);
		std::stringstream paramStr(paramsAsString);
		try {
		MSG (TrInfo, "2D copula (" << std::setw((int) floor(log10(nVars)) + 1) << i+1
		             << "," << std::setw((int) floor(log10(nVars)) + 1) << j+1
		             << ") of type '" << copType << "'");
		p2tgCop.reset(Copula2D::Cop2DInfo::make_new(copType, paramStr));
/*
		if (copNameMap.count(copType) == 0) {
			cerr << "Unknown 2D copula type `" << copType << "' .. aborting!" << endl;
			cout << "Known copula types are:";
			for (auto cIt = copNameMap.begin(); cIt != copNameMap.end(); ++ cIt)
				cout << " " << cIt->first;
			cout << endl;
			exit(1);
		}
		// copula type is OK -> it is safe to use copNameMap[copType]
		MSG (TrInfo, "2D copula (" << std::setw(floor(log10(nVars)) + 1) << i+1
		             << "," << std::setw(floor(log10(nVars)) + 1) << j+1
		             << ") of type '" << copType << "'");
		double dblP1, dblP2;
		unsigned uintP1;
		try {
		// NB: Cop2DTypeID = {indep, Clayton, Gumbel, Frank, Nelsen2, Nelsen18,
		//                    MarshallOlkin, data, normal, student};
		switch (copNameMap[copType]) {
		case Copula2D::Cop2DTypeID::indep: // independent margins
			p2tgCop = boost::make_shared<Copula2D::Cop2DIndep>();
			break;
		case Copula2D::Cop2DTypeID::Clayton: // independent margins
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DClayton>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::Gumbel: // independent margins
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DGumbel>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::Frank: // independent margins
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DFrank>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::Nelsen2: // independent margins
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DNelsen2>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::Nelsen18: // independent margins
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DNelsen18>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::MarshallOlkin: // independent margins
			tgFStr >> dblP1 >> dblP2; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DMarshallOlkin>(dblP1, dblP2);
			break;
		case Copula2D::Cop2DTypeID::data: // "sample"
			throw std::logic_error("2D sample copula is not yet implemented!");
			break;
		case Copula2D::Cop2DTypeID::normal: // "normal"
			tgFStr >> dblP1; // correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DNormal>(dblP1);
			break;
		case Copula2D::Cop2DTypeID::student: // student's t-copula
			tgFStr >> uintP1 >> dblP1; // DOF and correlation
			p2tgCop = boost::make_shared<Copula2D::Cop2DStudent>(uintP1, dblP1);
			break;
		default:
			cerr << "ERROR: file " << __FILE__ << ", line " << __LINE__
				  << " .. should never be here!" << endl;
			exit(1);
		}
*/
		}
		catch(std::exception& e) {
			cerr << "\nERROR in initialization of 2D copula no. " << n2Dcops + 1
			     << '\n' << "The error message was: " << e.what() << '\n';
			exit(1);
		}
		attach_2d_target(p2tgCop, i, j);
	}
}

