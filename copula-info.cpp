#include <iostream>
#include <fstream>

#include "common.hpp"
#include "copula-info.hpp"

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
	cMap["mixed"] = CopTypeID::mixed;
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

	throw std::logic_error("not yet implemented!");
}

