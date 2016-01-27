/** \file
	This is a driver illustrating the use of the cop-gen library.

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
	cout << "X = " << X << endl;
}

{
	cout << endl << "Example 2 - using sample as a target distribution" << endl;

	// the target sample - would normally be read from somewhere
	ublas::matrix<double> tgData(20,3); // variables in columns
	tgData <<=
		23.9, 17.5, 11.6,
		28.8, 22.3, 11.4,
		25.9, 18.5, 6.3,
		28.1, 20.3, 9.9,
		41.1, 27.5, 13.7,
		29.4, 12.7, 9.0,
		34.1, 21.9, 10.4,
		30.6, 21.0, 11.5,
		31.9, 21.3, 10.2,
		18.9, 16.9, 8.9,
		27.3, 20.5, 8.6,
		34.2, 26.2, 6.5,
		32.9, 25.5, 6.4,
		36.8, 29.5, 6.4,
		33.3, 25.9, 3.7,
		33.6, 25.7, 6.4,
		32.4, 27.2, 10.2,
		30.4, 23.1, 6.4,
		33.1, 25.7, 7.9,
		30.5, 24.5, 17.1;
	cout << "tgData = " << tgData << endl;

	uint nSc = 10;
	ublas::matrix<double> X(nSc, tgData.size2()); // variables in columns
	cout << "generating " << nSc << " scenarios .. "; std::cout.flush();
	double dist = ScenGenCop::gen_scen_sample(tgData, nSc, X);
	cout << "done." << endl;
	cout << "dist = " << dist << endl;
	cout << "X = " << X << endl;
}

	return 0;
}
