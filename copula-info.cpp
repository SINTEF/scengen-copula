#include <iostream>
#include <fstream>

#include "common.hpp"
#include "copula-info.hpp"

using namespace CopulaDef;

// ---------------------------------------------------------------------------
// class CopInfoBy2D


// ---------------------------------------------------------------------------
// class CopInfoNormal

CopInfoNormal::CopInfoNormal(UMatrix const & correls)
: CopInfoBy2D(correls.size1(), false), correlMat(correls)
{}

CopInfoNormal::CopInfoNormal(std::string const & tgFName)
: CopInfoBy2D(0, false)
{
	try {
		read_correl_mat(tgFName);
	}
	catch(std::exception& e) {
		std::cerr << e.what() << std::endl;
		throw; // re-throw the exception
	}
}


void CopInfoNormal::read_correl_mat(std::string const & tgFName)
{
	// read the input file
	std::ifstream tgCorrF(tgFName.c_str());
	if (!tgCorrF) {
		throw std::ios_base::failure("Could not open input file "
		                             + tgFName + "!");
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

