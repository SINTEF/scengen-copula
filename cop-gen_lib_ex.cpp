/// \file cop-gen_lib_ex.cpp
/**
	This is an example driver for the cop-gen library,
	illustrating how it can be used.

	Each example is put into its own scope, to make it stand-alone.

	\note The library needs the boost uBlas library for vectors and arrays
**/

#include "cop-gen_lib.hpp"
#include <boost/numeric/ublas/assignment.hpp> // assignment operator <<=
#include <boost/numeric/ublas/io.hpp>         // output operator <<

namespace ublas = boost::numeric::ublas; // shortcut name
using std::cout;
using std::endl;
using uint = unsigned; // shortcut name


int main(int argc, char *argv[]) {

{
	cout << "Example 1 - normal distribution" << endl;
	uint nVar = 3;
	ublas::vector<double> mean(nVar);
	ublas::vector<double> stD(nVar);
	ublas::matrix<double> corr(nVar, nVar); // alt: symmetric_matrix

	mean <<= 0, 1, 5;
	stD <<= 1, 1, 2;
	corr <<=  1 , 0.5, 0,
	         0.5,  1 , 0,
	         0.0, 0.0, 1;
	cout << "mean = " << mean << endl;
	cout << "st.D = " << stD << endl;
	cout << "corr = " << corr << endl;

	uint nSc = 5;
	ublas::matrix<double> X(nSc, nVar); // dim is not required
	cout << "generating " << nSc << " scenarios .. "; std::cout.flush();
	double dist = ScenGenCop::gen_scen_normal(mean, stD, corr, nSc, X, true, 1);
	cout << "done." << endl;
	cout << "dist = " << dist << endl;
	cout << "X = " << X;
}

	return 0;
}
