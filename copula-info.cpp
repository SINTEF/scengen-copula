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
				assert (p2Info2D(i, j).get() != NULL && "testing the null-check");
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
				assert (p2Info2D(i, j).get() != NULL && "testing the null-check");
				p2Info2D(i, j)->init_cdf_grid(gridPos);
			}
		}
	}
}
*/


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
		--i; --j; // assume the inputs count from 1, not from zero
		std::string copType;
		tgFStr >> copType;
		if (! tgFStr.good())
			continue;

		// NEW - NOT YET FINISHED!
		std::string paramsAsString;
		std::getline(tgFStr, paramsAsString);
		std::stringstream paramStr(paramsAsString);
		try {
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
			cerr << "Error: There was some problem initializing bivariate copulas!"
			     << endl << "       The error message was: " << e.what() << endl;
			exit(1);
		}
		attach_2d_target(p2tgCop, i, j);
	}
}

