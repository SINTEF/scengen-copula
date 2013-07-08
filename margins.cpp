#include <fstream>
#include <ios>

#include "margins.hpp"

using namespace MarginDistrib;

// ---------------------------------------------------------------------------
// generic routines

void MarginDistrib::make_marg_name_map(MargNameMapT & mMap) {
	// note: when listing the map, it goes from last to first
	mMap["sample"] = MargTypeID::sample;
	mMap["normal"] = MargTypeID::normal;
	mMap["fixed"] = MargTypeID::fixed;
	mMap["mixed"] = MargTypeID::mixed;
}


// ----------------------------------------------------------------------------
// class MarginsInfo

void MarginsInfo::get_margin_distr(MatrixI const & ranks, MatrixD & values)
{
	DimT i;//, j;
	assert (ranks.size1() == nM && "sanity check");
	DimT nVals = ranks.size2();
	// We cannot send the row as the target -> send a vector & copy the results.
	// Note that we do not have to clear the vector, as it gets overwritten..
	VectorD resVect(nVals);
	values.resize(ranks.size1(), ranks.size2());
	for (i = 0; i < nM; ++i) {
		p2margins(i)->inv_cdf(ublas::row(ranks, i), resVect);
		assert (resVect.size() == nVals && "sanity check");
		ublas::row(values, i) = resVect;
		#ifndef NDEBUG
			for (DimT j = 0; j < nVals; ++j) {
				assert(isEq(values(i, j), resVect(j)) && "checking row assignment");
			}
		#endif
	}
}


// ----------------------------------------------------------------------------
// class MixedMargins

// attach margin
void MixedMargins::attach_margin(UnivarMargin::Ptr & p2marg, DimT const index)
{
	p2margins[index] = p2marg;
}


// ----------------------------------------------------------------------------
// class NormalMargins

NormalMargins::NormalMargins(VectorD const & EVs, VectorD const & SDs,
                             UnivarMargin::SamplePP const postP)
: MarginsInfo(EVs.size()), means(EVs), stDevs(SDs)
{
	init_margins(postP);
}

NormalMargins::NormalMargins(std::string fName,
                             UnivarMargin::SamplePP const postP)
: MarginsInfo(0)
{
	read_from_file(fName);
	init_margins(postP);
}

void NormalMargins::read_from_file(std::string fName)
{
	// read the input file
	std::ifstream inFile(fName.c_str());
	if (!inFile) {
		throw std::ios_base::failure("Could not open input file " + fName + "!");
	}
	// the first datum is the number of variables
	inFile >> nM;
	if (nM <= 0) {
		throw std::range_error("number of margins must be positive");
	}
	p2margins.resize(nM);
	means.resize(nM);
	stDevs.resize(nM);

	// the rest of the file consists of means and standard deviations,
	// with one margin per line
	DimT i;
	for (i = 0; i < nM; ++i) {
		inFile >> means(i) >> stDevs(i);
		if (stDevs(i) < DblEps) {
			throw std::range_error("standard dev. must be positive");
		}
	}
}


void NormalMargins::init_margins(UnivarMargin::SamplePP const postP)
{
	DimT i;
	for (i = 0; i < nM; ++i) {
		p2margins(i).reset(new MarginNormal(means(i), stDevs(i), postP));
	}

}


// ----------------------------------------------------------------------------
// class SampleMargins

SampleMargins::SampleMargins(MatrixD const & samples,
                             UnivarMargin::SamplePP const postP)
: MarginsInfo(samples.size1()), data(samples)
{
	init_margins(postP);
}


void SampleMargins::init_margins(UnivarMargin::SamplePP const postP)
{
	DimT i;
	for (i = 0; i < nM; ++i) {
		p2margins(i).reset(new MarginSample(ublas::row(data, i), postP));
	}
}
