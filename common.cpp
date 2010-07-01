#include "common.hpp"
#include "ranker.h"

void get_ranks(const std::vector<double> &inputVect, std::vector<int> &ranks) {
	// using the rank() from "ranker.h"
	// this gives numbers from 1 to N, so we have to subtract 1!
	rank(inputVect, ranks);
	unsigned len = ranks.size();
	for (unsigned i = 0; i < len; i++) {
		ranks[i]--;
	}
}

/*
void get_ranks(double const inputVect[], std::vector<int> ranks); {
	size_t len = sizeof(inputVect) / sizeof(double);
	rank(inputVect, len, ranks);
	for (unsigned i = 0; i < len; i++) {
		ranks[i]--;
	}
}
*/
