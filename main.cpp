#include <iostream>
//#include <ctime>

#include "cop2Dsample.hpp"
#include "copula-sample.hpp"

using namespace std;
using namespace Copula2D;
using namespace CopulaScen;

int main(int argc, char *argv[]) {
	// set the random seed
  //srand ( time(NULL) );
  srand(0); // for debugging

	const int defaultNmbSamples = 10;
	int nS = (argc > 1 ? atoi(argv[1]) : defaultNmbSamples);

	Cop2DClayton tgCopCl(-0.75);
	Cop2DIndep tgCopInd;
	Cop2DNelsen18 tgCopN18(5.0);
	Cop2DMarshallOlkin tgCopMO(0.5, 0.75);
	Cop2DInfTr tgCopTr(&tgCopMO);

	// using the bivariate code directly
	/*
		NB: This gives different results, since this codes generates the copula
		    by columns, while the general one does it by columns!
	*/
	Cop2DSample bivarCopSc(nS, &tgCopMO);
	bivarCopSc.gen_heur();
	bivarCopSc.print_as_txt("test.txt");


	// using the new multivariate code
	/*
	CopulaSample copSc(2, nS);
	copSc.attach_tg_2Dcop(&tgCopTr, 0, 1);
	copSc.gen_sample();
	copSc.print_as_txt("out_cop.txt");
	copSc.print_2D_as_txt(0, 1, "out_2D-cop_0-1.txt");
	*/


	// using the new multivariate code - 3D
	CopulaSample copSc(3, nS);
	copSc.attach_tg_2Dcop(&tgCopTr, 0, 1);
	//copSc.attach_tg_2Dcop(&tgCopInd, 0, 2);
	copSc.attach_tg_2Dcop(&tgCopMO, 0, 2);
	//copSc.attach_tg_2Dcop(&tgCopN18, 1, 2);
	copSc.attach_tg_2Dcop(&tgCopInd, 1, 2);
	copSc.gen_sample();
	copSc.print_as_txt("out_cop.txt");
	copSc.print_2D_as_txt(0, 1, "out_2D-cop_0-1.txt");
	copSc.print_2D_as_txt(0, 2, "out_2D-cop_0-2.txt");
	copSc.print_2D_as_txt(1, 2, "out_2D-cop_1-2.txt");


	return 0;
}
