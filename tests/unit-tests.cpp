#include "../common.hpp"

#include <gtest/gtest.h>
#include <iostream>

OutputLevel outLvl = TrError;

GTEST_API_ int main(int argc, char **argv) {
	std::cout << "Running unit tests for copula scen-gen" << std::endl;

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
