#include <fstream>
#include <ios>
#include <map>

#include "margins.hpp"
#include "margin-distrib_moments.hpp"

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
	mMap["moments"] = MargTypeID::moments;
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
MixedMargins::MixedMargins(std::string const & tgFName, DimT const nScens,
                           DimT const nVars)
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
		MSG (TrInfo, "margin " << std::setw(floor(log10(nM)) + 1) << i+1
		             << " of type '" << distrID << "'");

		// distribution type is OK -> it is safe to use distrNameMap[distrID]
		try {
			UnivarMargin::Ptr p2marg;
			switch (distrNameMap[distrID]) {
			case MargDistribID::moments: // "moments"
				p2marg = boost::make_shared<MarginMoments>(paramStr, nScens);
				break;
			case MargDistribID::normal: // "normal"
				p2marg = boost::make_shared<MarginNormal>(paramStr);
				break;
			case MargDistribID::sample: // sample-based
				p2marg = boost::make_shared<MarginSample>(paramStr, nScens);
				break;
			case MargDistribID::triang: // triangular
				p2marg = boost::make_shared<MarginTriang>(paramStr, false);
				break;
			case MargDistribID::triangX: // triangular with values at min & max
				p2marg = boost::make_shared<MarginTriang>(paramStr, true);
				break;
			case MargDistribID::exponential: // exponential
				p2marg = boost::make_shared<MarginExp>(paramStr);
				break;
			case MargDistribID::beta: // exponential
					p2marg = boost::make_shared<MarginBeta>(paramStr);
				break;
			case MargDistribID::lognormal: // exponential
					p2marg = boost::make_shared<MarginLognormal>(paramStr);
				break;
			case MargDistribID::poisson: // exponential
					p2marg = boost::make_shared<MarginPoisson>(paramStr);
				break;
			case MargDistribID::uniform: // continuous uniform distribution
					p2marg = boost::make_shared<MarginUniformC>(paramStr);
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


// --------------------------------------------------------------------------
// class MarginsByMoments

MarginsByMoments::MarginsByMoments(std::string fName, DimT const nScen,
                                   int const FoM, bool const matFmt)
: MarginsInfo(0), formOfMoments(FoM), nSc(nScen)
{
	read_from_file(fName, matFmt);
	init_margins();
}

void MarginsByMoments::read_from_file(std::string fName, bool const matFmt)
{
	// read the input file
	std::ifstream inFile(fName.c_str());
	if (!inFile) {
		throw std::ios_base::failure("Could not open input file " + fName + "!");
	}
	if (matFmt) {
		// using the (old) matrix-style input format:
		// 	4	dim
		// 	one moment per line (4 lines)
		inFile >> nM;
		if (nM != 4)
			throw std::runtime_error("the matrix-style format for moments implies "
			                         "that the first number is 4");
	} else {
		// the new format is:
		// 	dim
		// 	one margin per line (4 columns)
	}
	inFile >> nM;
	if (nM <= 0) {
		throw std::range_error("number of margins must be positive");
	}
	p2margins.resize(nM);
	tgMoments.resize(nM);

	// read the values
	DimT i, j;
	if (matFmt) {
		// one moment per line (4 lines)
		for (i = 0; i < nM; ++i)
			tgMoments[i].resize(4);
		for (j = 0; j < 4; ++j) {
			for (i = 0; i < nM; ++i) {
				inFile >> tgMoments[i](j);
			}
		}
	} else {
		// one margin per line (4 columns)
		for (i = 0; i < nM; ++i) {
			tgMoments[i].resize(4);
			for (j = 0; j < 4; ++j)
				inFile >> tgMoments[i](j);
		}
	}
}

void MarginsByMoments::init_margins()
{
	DimT i;
	for (i = 0; i < nM; ++i) {
		p2margins(i).reset(new MarginMoments(tgMoments[i], nSc, formOfMoments));
	}

}
