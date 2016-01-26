//#include "../cop-gen_lib.hpp"
#include "../copula-sample.hpp"
#include "../margins.hpp"

#include <boost/numeric/ublas/assignment.hpp>
#include <gtest/gtest.h>

using namespace Copula2D;

/// simple example of a normal-distribution target
class NormalDistribEx : public ::testing::Test {
public:
	NormalDistribEx( ) {
		// initialization code here
		nVar = 4;
		tgMean.resize(nVar);
		tgStD.resize(nVar);
		tgCorr.resize(nVar, nVar);

		tgMean <<= 1.0, 2.0, 3.0, 4.0;
		tgStD  <<= 1.0, 1.0, 2.0, 2.0;
		tgCorr <<= 1.0, 0.5, 0.3, 0.3,
				   0.5, 1.0, 0.3, 0.3,
				   0.3, 0.3, 1.0, 0.5,
				   0.3, 0.3, 0.5, 1.0;
		nSc = 10;
	}

	void SetUp( ) {
		// code here will execute just before the test ensues
		p2tgCop = boost::make_shared<CopulaDef::CopInfoNormal>(tgCorr);
		p2tgMargins
			= boost::make_shared<MarginDistrib::NormalMargins>(tgMean, tgStD);
	}

	void TearDown( ) {
		// code here will be called just after the test completes
		// ok to through exceptions from here if need be
	}

	~NormalDistribEx()  {
		// cleanup any pending stuff, but no exceptions allowed
	}

	// put in any custom data members that you need
	VectorD tgMean, tgStD;
	MatrixD tgCorr;
	DimT nVar, nSc;
	CopulaDef::CopInfoNormal::Ptr p2tgCop;
	MarginDistrib::NormalMargins::Ptr p2tgMargins;
};


/// test generating scenarios from normal distribution
TEST_F(NormalDistribEx, TestCopulaGen) {
	// generate scenarios for the copula (as ranks)
	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	double scCopDist = copSc.gen_sample();
	EXPECT_NEAR(scCopDist, 0.0, 1e-10)
	            << "the copula generation error was " << scCopDist;
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);
	EXPECT_EQ(copRanks.size1(), nVar) << "checking the number of variables";
	EXPECT_EQ(copRanks.size2(), nSc) << "checking the number of scenarios";

	// transform margins to the target distribution
	MatrixD scens;
	p2tgMargins->get_margin_distr(copRanks, scens);
	EXPECT_EQ(scens.size1(), nVar) << "checking the number of variables";
	EXPECT_EQ(scens.size2(), nSc) << "checking the number of scenarios";
}


/// test fixing the first margin
TEST_F(NormalDistribEx, TestFixedMargins1) {
	// generate scenarios for the copula (as ranks)
	CopulaScen::CopulaSample origCopSc(p2tgCop, nSc);
	origCopSc.gen_sample();
	MatrixI origCopRanks;
	origCopSc.get_result_ranks(origCopRanks);

	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	MatrixI fixedMargs(1, nSc);
	ublas::row(fixedMargs, 0) = ublas::row(origCopRanks, 0);
	copSc.fix_marg_values(fixedMargs);

	double scCopDist = copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);

	for (unsigned s = 0; s < nSc; ++s)
		EXPECT_EQ(copRanks(0, s), origCopRanks(0, s)) << "check fixed values";

	EXPECT_NEAR(scCopDist, 0.0, 1e-6)
				<< "the copula generation error was " << scCopDist;
}

/// test fixing the first margin with values from the second (i.e., not sorted)
TEST_F(NormalDistribEx, TestFixedMargins2) {
	// generate scenarios for the copula (as ranks)
	CopulaScen::CopulaSample origCopSc(p2tgCop, nSc);
	origCopSc.gen_sample();
	MatrixI origCopRanks;
	origCopSc.get_result_ranks(origCopRanks);

	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	MatrixI fixedMargs(1, nSc);
	ublas::row(fixedMargs, 0) = ublas::row(origCopRanks, 1);
	copSc.fix_marg_values(fixedMargs);

	double scCopDist = copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);

	for (unsigned s = 0; s < nSc; ++s)
		EXPECT_EQ(copRanks(0, s), origCopRanks(1, s)) << "check fixed values";

	EXPECT_NEAR(scCopDist, 0.0, 1e-6)
				<< "the copula generation error was " << scCopDist;
}

/// test fixing the first two margins, reversed
TEST_F(NormalDistribEx, TestFixedMargins3) {
	// generate scenarios for the copula (as ranks)
	CopulaScen::CopulaSample origCopSc(p2tgCop, nSc);
	origCopSc.gen_sample();
	MatrixI origCopRanks;
	origCopSc.get_result_ranks(origCopRanks);

	CopulaScen::CopulaSample copSc(p2tgCop, nSc);
	MatrixI fixedMargs(2, nSc);
	ublas::row(fixedMargs, 0) = ublas::row(origCopRanks, 1);
	ublas::row(fixedMargs, 1) = ublas::row(origCopRanks, 0);
	copSc.fix_marg_values(fixedMargs);

	//outLvl = TrDetail;

	double scCopDist = copSc.gen_sample();
	MatrixI copRanks;
	copSc.get_result_ranks(copRanks);

	for (unsigned s = 0; s < nSc; ++s) {
		EXPECT_EQ(copRanks(0, s), origCopRanks(1, s)) << "check fixed values";
		EXPECT_EQ(copRanks(1, s), origCopRanks(0, s)) << "check fixed values";
	}

	EXPECT_NEAR(scCopDist, 0.0, 1e-6)
				<< "the copula generation error was " << scCopDist;
}
