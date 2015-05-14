#include "../common.hpp"
#include "../cop2Dinfo.hpp"
#include <gtest/gtest.h>

using namespace Copula2D;

OutputLevel outLvl = TrInfo;

TEST(Cop2DGridTest, CdfR) {
	int const nSc = 5;

	Cop2DInfo::Ptr p2C = boost::make_shared<Cop2DIndep>();
	p2C->set_nmb_scens(nSc);
	Cop2DInfo::Ptr p2G = boost::make_shared<Cop2DGrid>(p2C);

	int i, j;
	for (i = 0; i < nSc; ++i) {
		for (j = 0; j < nSc; ++j) {
			EXPECT_DOUBLE_EQ(p2C->cdfR(i,j), p2G->cdfR(i,j))
				<< "cdfR differs for (i,j) = (" << i << "," << j << ")";
		}
	}
}


TEST(Cop2DApproxGridTest, ExactCdfR) {
	int const nSc = 10;
	int const gridN = nSc / 2;

	Cop2DInfo::Ptr p2C = boost::make_shared<Cop2DIndep>();
	p2C->set_nmb_scens(nSc);

	Cop2DInfo::Ptr p2G = boost::make_shared<Cop2DGrid>(p2C, gridN);
	int i, j;
	for (i = 1; i < nSc; i = i + 2) {
		for (j = 1; j < nSc; j = j + 2) {
			EXPECT_DOUBLE_EQ(p2C->cdfR(i,j), p2G->cdfR(i,j))
				<< "cdfR differs for (i,j) = (" << i << "," << j << ")";
		}
	}
}


// This uses the fact that the interpolation is exact for independent copula!
TEST(Cop2DApproxGridTest, InterpolatedCdfR) {
	int const nSc = 10;
	int const gridN = nSc / 2;

	Cop2DInfo::Ptr p2C = boost::make_shared<Cop2DIndep>();
	p2C->set_nmb_scens(nSc);

	Cop2DInfo::Ptr p2G = boost::make_shared<Cop2DGrid>(p2C, gridN);
	int i, j;
	for (i = 0; i < nSc; i = i + 2) {
		for (j = 0; j < nSc; j = j + 2) {
			EXPECT_NEAR(p2C->cdfR(i,j), p2G->cdfR(i,j), 1e-9)
				<< "cdfR differs for (i,j) = (" << i << "," << j << ")";
		}
	}
}
