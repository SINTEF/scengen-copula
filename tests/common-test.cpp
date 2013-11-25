#include "../common.hpp"
#include <gtest/gtest.h>

//! this file uses C++11 features (vector assignment)!
//! - with gcc, using -std=c++0x or -std=c++11 causes erros in gtest-port.h
//! - as a workaround, one can use -std=gnu++0x or -std=gnu++11 instead,
//!   see http://stackoverflow.com/questions/6312151/


/// tests for objects defined in common.hpp
TEST(TestCommon, RankFunctions) {
	std::vector<double> valV;
	std::vector<int> rankV;
	std::vector<int> expRV;

	valV = {1.5, 0.5, 1.0, -1};
	expRV = {3, 1, 2, 0};
	get_ranks(valV, rankV);
	EXPECT_EQ (expRV, rankV) << "get_ranks() failed a basic test";

	valV = {1.5, 0.5, 1.0, 0.5, 0, 0.5};
	expRV = {5, 1, 4, 2, 0, 3}; // ties should be resolved based on order
	get_ranks(valV, rankV);
	EXPECT_EQ (expRV, rankV) << "get_ranks() failed a tie-resolving test";

	/// \todo test get_ranks_or_rows()

	EXPECT_EQ(0.0, rank2U01(-1, 13)) << "rank2U01() failed test with r=-1";
	EXPECT_EQ(1.0, rank2U01(12, 13)) << "rank2U01() failed test with r=N-1";

	EXPECT_EQ( -1, u012Rank(0.0, 13)) << "u012Rank() failed test with z=0";
	EXPECT_EQ( 12, u012Rank(1.0, 13)) << "u012Rank() failed test with z=1";

	EXPECT_EQ(6, u012Rank(rank2U01(6, 13), 13))
		<< "failed u012Rank(rank2U01()) inverse test with r=-1";
	EXPECT_EQ(-1, u012Rank(rank2U01(-1, 13), 13))
		<< "failed u012Rank(rank2U01()) inverse test with r=-1";
	EXPECT_EQ(12, u012Rank(rank2U01(12, 13), 13))
		<< "failed u012Rank(rank2U01()) inverse test with r=-1";
	EXPECT_EQ(0.0, rank2U01(u012Rank(0.0, 13), 13))
		<< "failed rank2U01(u012Rank()) inverse test with z=0";
	EXPECT_EQ(0.5, rank2U01(u012Rank(0.5, 14), 14)) // 0.5 needs even N
		<< "failed rank2U01(u012Rank()) inverse test with z=0.5";
	EXPECT_EQ(1.0, rank2U01(u012Rank(1.0, 13), 13))
		<< "failed rank2U01(u012Rank()) inverse test with z=1";

	/// \todo test fix_mean_std()
}
