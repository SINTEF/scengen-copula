#include <fstream>
#include <ios>

#include "margins.hpp"

using namespace MarginDistrib;
using std::cout;
using std::endl;
using std::cerr;

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

// constructor from a file
MixedMargins::MixedMargins(std::string const & tgFName, DimT const nVars)
: MarginsInfo(nVars)
{
	std::ifstream margSpecF(tgFName);
	if (!margSpecF) {
		throw std::ios_base::failure("Could not open input file `" + tgFName
		                             + "'!");
	}
	margSpecF >> nM;
	if (nVars > 0 && nVars != nM) {
		std::stringstream errMsgStr;
		errMsgStr << "Wrong number of margins in file `" << tgFName << "' "
		          << "- expected " << nVars << ", got " << nM << "!";
		throw std::runtime_error(errMsgStr.str());
	}
	p2margins.resize(nM);

	if (distrNameMap.size() == 0) {
		make_distrib_name_map(distrNameMap);
	}

	std::string distrID;
	std::string paramsAsString;
	for (DimT i =0; i < nM; ++i) {
		margSpecF >> distrID;
		std::getline(margSpecF, paramsAsString);
		std::stringstream paramStr(paramsAsString);

		if (distrNameMap.count(distrID) == 0) {
			cout << "Known distribution types are:";
			for (auto dIt = distrNameMap.begin(); dIt != distrNameMap.end(); ++ dIt)
				cout << " " << dIt->first;
			cout << endl;
			throw std::runtime_error("Unknown marginal distribution type `"
			                         + distrID + "'");
		}

		// distribution type is OK -> it is safe to use distrNameMap[distrID]
		try {
			UnivarMargin::Ptr p2marg;
			switch (distrNameMap[distrID]) {
			case MargDistribID::moments: // "sample"
				MSG (TrInfo2, "margin of type 'sample'");
				throw std::logic_error("not yet implemented");
				break;
			case MargDistribID::normal: // "normal"
				MSG (TrInfo, "margin of type 'normal'");
				throw std::logic_error("not yet implemented");
				break;
			case MargDistribID::triang: // generic mixed of 2D copulas
				MSG (TrInfo, "margin of type 'triangular'");
				// create a new object of the specific class
				p2marg = boost::make_shared<MarginTriang>(paramStr, false);
				break;
			case MargDistribID::triangX: // generic mixed of 2D copulas
				MSG (TrInfo, "margins of type 'triangular with extremes'");
				// create a new object of the specific class
				p2marg = boost::make_shared<MarginTriang>(paramStr, true);
				break;
			default:
				throw std::logic_error("Should never be here!");
			}
			attach_margin(p2marg, i);
		}
		catch(std::exception& e) {
			cerr << "Error: There was some problem initializing a margin!" << endl
				  << "       The error message was: " << e.what() << endl;
			exit(1);
		}
	}
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
