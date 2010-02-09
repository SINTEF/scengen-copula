#include <iostream>
//#include <ctime>

#include "cop2Dsample.hpp"
#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;

int main(int argc, char *argv[]) {
	// set the random seed
  srand ( time(NULL) );
  //srand(0); // for debugging

	const int defaultNmbSamples = 10;
	int nS = (argc > 1 ? atoi(argv[1]) : defaultNmbSamples);

	//Cop2DClayton tgCop(-0.75);
	//Cop2DIndep tgCop;
	//Cop2DNelsen18 tgCop(5.0);
	Cop2DMarshallOlkin tgCop(0.5, 0.75);

	// using the bivariate code directly
	Cop2DSample bivarCopSc(nS, &tgCop);
	bivarCopSc.gen_heur();
	bivarCopSc.print_as_txt("test.txt");

	// using the new multivariate code
	CopulaSample copSc(2, nS);
	copSc.attach_tg_2Dcop(&tgCop, 0, 1);
	copSc.gen_sample();
	copSc.print_as_txt("out_cop.txt");
	copSc.print_2D_as_txt(0, 1, "out_2D-cop_0-1.txt");

	return 0;
}
