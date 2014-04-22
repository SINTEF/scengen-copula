// implementation of class CopInfoNormal

#include <iostream>
#include <fstream>

#include "copula-info.hpp"

using namespace std;
using namespace CopulaDef;


CopInfoNormal::CopInfoNormal(MatrixD const & correls)
: CopInfoBy2D(correls.size1(), false), correlMat(correls)
{
	setup_2d_targets(); // create the matrix of bivariate copula objects
}

CopInfoNormal::CopInfoNormal(std::string const & tgFName)
: CopInfoBy2D(0, false)
{
	try {
		read_correl_mat(tgFName);
	}
	catch(std::exception& e) {
		cerr << "Error: while reading data file `" << tgFName << "'!" << endl;
		//cerr << "       The error message was: " << e.what() << endl;
		throw; // re-throw the exception
	}
	setup_2d_targets(); // create the matrix of bivariate copula objects
}


void CopInfoNormal::read_correl_mat(std::string const & tgFName)
{
	// read the input file
	std::ifstream tgCorrF(tgFName.c_str());
	if (!tgCorrF) {
		throw std::ios_base::failure("Could not open input file `"
		                             + tgFName + "'!");
	}
	get_correl_matrix_from_stream(tgCorrF, correlMat);

	tgCorrF.close();
}


void CopInfoNormal::setup_2d_targets()
{
	assert (nVars > 0 && "must have a known dimension by now");
	if (p2Info2D.size1() != nVars || p2Info2D.size2() != nVars) {
		if (p2Info2D.size1() * p2Info2D.size2() > 0) {
			cerr << "Warning: resizing the matrix of 2D-copula objects!" << endl;
		}
		p2Info2D.resize(nVars, nVars);
	}

	DimT i, j;
	for (i = 0; i < nVars; i++) {
		for (j = i+1; j < nVars; j++) {
			attach_2d_target(new Copula2D::Cop2DNormal(correlMat(i, j)), i, j);
		}
	}
}
