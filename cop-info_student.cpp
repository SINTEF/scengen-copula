/** \file
	implementation of class CopInfoStudent
**/

#include "copula-info.hpp"

#include <iostream>
#include <fstream>

using namespace std;
using namespace CopulaDef;


CopInfoStudent::CopInfoStudent(unsigned degF, MatrixD const & correls)
: CopInfoNormal(correls), dof(degF)
{
	setup_2d_targets(); // create the matrix of bivariate copula objects
}

CopInfoStudent::CopInfoStudent(unsigned degF, std::string const & tgFName)
: CopInfoNormal(tgFName), dof(degF)
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

CopInfoStudent::CopInfoStudent(std::string const & tgFName)
: CopInfoNormal(tgFName), dof(0)
{
	try {
		get_params_from_file(tgFName);
	}
	catch(std::exception& e) {
		cerr << "Error: while reading data file `" << tgFName << "'!" << endl;
		//cerr << "       The error message was: " << e.what() << endl;
		throw; // re-throw the exception
	}
	setup_2d_targets(); // create the matrix of bivariate copula objects
}



void CopInfoStudent::setup_2d_targets()
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
			attach_2d_target(boost::make_shared<Copula2D::Cop2DStudent>(dof, correlMat(i, j)), i, j);
		}
	}
}


void CopInfoStudent::get_params_from_file(std::string const & tgFName)
{
	// read the input file
	std::ifstream paramF(tgFName.c_str());
	if (!paramF) {
		throw std::ios_base::failure("Could not open input file `"
		                             + tgFName + "'!");
	}
	paramF >> dof; // first number is the dof
	get_correl_matrix_from_stream(paramF, correlMat);

	paramF.close();
}
