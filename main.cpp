#include <iostream>
//#include <ctime>

#include "cop2Dsample.hpp"

using namespace std;
using namespace Copula2D;

int main(int argc, char *argv[]) {
	// set the random seed
  srand ( time(NULL) );
  //srand(0); // for debugging

	const int defaultNmbSamples = 10;
	int nS = (argc > 1 ? atoi(argv[1]) : defaultNmbSamples);

	//Cop2DClayton tgCop(-0.75);
	//Cop2DIndep tgCop;
	Cop2DNelsen18 tgCop(5.0);
	Cop2DSample  copSc(nS, &tgCop);

	copSc.gen_heur();
	copSc.print_as_txt("test.txt");

	return 0;
}
