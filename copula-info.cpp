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
	cMap["sample"] = cSample;
	cMap["normal"] = cNormal;
	cMap["indep"] = cIndep;
}


// ---------------------------------------------------------------------------
// class CopIndep

double CopIndep::cdf(VectorD const u) const {
	assert (u.size() == nVars && "dimension check");

	double F = u[0];
	DimT i;
	for (i = 1; i < nVars; ++i)
		F *= u[i];

	return F;
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

