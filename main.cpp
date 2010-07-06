#include <iostream>
#include <fstream> // TEMP!
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

/*
  Cop2DSample::CandList cList(2, 1e-4);
	int lIdx; double lVal;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items" << endl;
	//
	lIdx = 1; lVal = 0.2;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;
	//
	lIdx = 2; lVal = 0.1;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;
	//
	lIdx = 3; lVal = 0.0;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;
	//
	lIdx = 4; lVal = 7.0;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;
	//
	lIdx = 5; lVal = 6.0;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;

	//
	lIdx = 13; lVal = 2.0;
	cList.insert_cand(lIdx, lVal);
	cout << "CandList test: added cand (" << lIdx << ", " << lVal << ")" << endl;
  cout << "CandList test: list includes " << cList.get_num_cand()
       << " items; the best one is (" << cList.get_best_cand_index()
       << ", " << cList.get_best_cand_value() << ")" << endl;
	return 0;
*/

	const int defaultNmbSamples = 10;
	int nS = (argc > 1 ? atoi(argv[1]) : defaultNmbSamples);

	// ------------- TEST -----------------------------
	ifstream testFile("Dorey_cop_2D.txt");
	if (!testFile) {
		cerr << "Problem with the input file!" << endl;
		exit(1);
	}
	int nRows, nCols;
	testFile >> nRows >> nCols;
	std::vector<double> margin[2];
	int i, j;
	for (i = 0; i < 2; i++)
		margin[i].resize(nRows);
	for (i = 0; i < nRows; i++)
		for (j = 0; j < nCols; j++)
			testFile >> margin[j][i];
	Cop2DData< std::vector<double> > tgCopData(&margin[0], &margin[1], nRows, nS);
	CopulaSample testCopSc(2, nS);
	testCopSc.attach_tg_2Dcop(&tgCopData, 0, 1);
	testCopSc.gen_sample();
	testCopSc.print_as_txt("test_out_cop.txt");
	cout << "target Dorey_cop_2D.txt - finished." << endl;
	//exit(0);
	// -------------------------------------------------

	Cop2DClayton tgCopCl(-0.75);
	Cop2DIndep tgCopInd;
	Cop2DNelsen18 tgCopN18(5.0);
	Cop2DMarshallOlkin tgCopMO(0.5, 0.75);
	Cop2DInfTr tgCopTr(&tgCopMO);

	// using the bivariate code directly
	/*
		NB: This gives different results, since this codes generates the copula
		    by columns, while the general one does it by columns!
	*//*
	Cop2DSample bivarCopSc(nS, &tgCopMO);
	bivarCopSc.gen_heur();
	bivarCopSc.print_as_txt("test.txt");
	*/

	// using the new multivariate code

	CopulaSample copSc(2, nS);
	copSc.attach_tg_2Dcop(&tgCopTr, 0, 1);
	copSc.gen_sample();
	copSc.print_as_txt("out_cop.txt");
	copSc.print_2D_as_txt(0, 1, "out_2D-cop_0-1.txt");


	/*
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
	*/

	return 0;
}
