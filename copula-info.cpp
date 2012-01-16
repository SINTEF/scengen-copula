#include <iostream>
#include <fstream>

#include "common.hpp"
#include "copula-info.hpp"

using namespace CopulaDef;
using std::cout;
using std::cerr;
using std::endl;

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
// class CopInfoNormal

CopInfoNormal::CopInfoNormal(MatrixD const & correls)
: CopInfoBy2D(correls.size1(), false), correlMat(correls)
{}

CopInfoNormal::CopInfoNormal(std::string const & tgFName)
: CopInfoBy2D(0, false)
{
	try {
		read_correl_mat(tgFName);
	}
	catch(std::exception& e) {
		cerr << "Error: Could not open data file `" << tgFName
		     << "' for the normal copula!" << endl;
		cerr << "       The error message was: " << e.what() << endl;
		throw; // re-throw the exception
	}
}


void CopInfoNormal::read_correl_mat(std::string const & tgFName)
{
	// read the input file
	std::ifstream tgCorrF(tgFName.c_str());
	if (!tgCorrF) {
		throw std::ios_base::failure("Could not open input file `"
		                             + tgFName + "'!");
	}
	tgCorrF >> correlMat; // this should read the matrix, including dimensions

	nVars = correlMat.size1();
	assert (correlMat.size2() == nVars && "must be a square matrix");

	#ifndef NDEBUG
		for (DimT i = 0; i < nVars; ++i) {
			if (! isEq(correlMat(i, i), 1.0)) {
				throw std::range_error("correl(i,i) must be 1");
			}
			for (DimT j = 0; j < i; ++j) {
				double rho = correlMat(i, j);
				if (rho < -1 || rho > 1 || !isEq(rho, correlMat(j, i))) {
					throw std::range_error("error in the correlation matrix");
				}
			}
		}
	#endif
}


void CopInfoNormal::setup_2d_targets()
{
	DimT i, j;
	p2Info2D.resize(nVars, nVars);
	for (i = 0; i < nVars; i++) {
		for (j = i+1; j < nVars; j++) {
			attach_2d_target(new Copula2D::Cop2DNormal(correlMat(i, j)), i, j);
		}
	}
}

