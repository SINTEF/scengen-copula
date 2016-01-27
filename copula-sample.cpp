// implementation of class CopulaSample

#include "copula-sample.hpp"

#include <iostream>
#include <fstream>

using namespace std;
using namespace CopulaScen;


// --------------------------------------------------------------------------
// CONSTRUCTORS AND DESTRUCTORS
CopulaSample::CopulaSample(CopulaDef::CopInfoBy2D::Ptr p2tg, DimT const S,
                           unsigned const numCandPts)
: nVar(p2tg->dim()), nSc(S),
  haveSc4Marg(nVar, false),
  p2copInfo(p2tg.get()),
  p2tgInfo(p2tg->get_pts_to_2d_targets()),
  p2sample2D(nVar, nVar),
  p2prob(nullptr), sample(nVar),
  minNumCandScens(numCandPts)
{
	DimT i, j;
	for (i = 0; i < nVar; i++) {
		for (j = 0; j < nVar; j++) {
			p2sample2D(i,j) = nullptr;
		}
		sample[i].resize(nSc);
	}

	//p2tg->init_cdf_grids(S); // done in the 2D classes

	// at the moment, only the transposed targets are allocated in the class
	//allocTgCops.reserve(nVar * (nVar-1) / 2);
}


double CopulaSample::gen_new_margin(DimT const marg)
{
	DimT i, j, s;
	DimT tg, iR;

	if (haveSc4Marg[marg])
		throw std::logic_error("trying to generate an existing margin");

	/// vector of target specifications for all the involved 2D copulas
	std::vector<Cop2DInfo const*> tg2Dcopulas;
	tg2Dcopulas.reserve(nVar);

	/// list of already existing margins that are being matched
	std::vector<int> oldMargins;
	oldMargins.reserve(nVar);

	double CdfDistEps = DblEps; // DblEps should give 'optimal' discretization

	// find the copulas to be matched
	// for the moment, assume that the new variable is the second index!!!
	for (i = 0; i < nVar; i++) {
		if (i != marg && haveSc4Marg[i] && p2tgInfo(i, marg)) {
			// store the margin in the list of copulas (i, marg) to match
			oldMargins.push_back(i);
			tg2Dcopulas.push_back(p2tgInfo(i, marg).get());
			// create a new copula-2D-sample object
			stringstream cop2DId; // using stringstream to get simple conversions
			cop2DId << "sample_" << i << "_" << marg;
			p2sample2D(i, marg) = new Cop2DSample(nSc, p2tgInfo(i, marg).get(),
			                                      cop2DId.str());
			TRACE (TrDetail2, "Created new Cop2DSample with id = "
			       << cop2DId.str());

			bool canHaveDuplicates = i < nmbMargsWithDuplicates;

			// initialize the 2D-sample generator with scenarios of the known margin
			p2sample2D(i, marg)->set_scen_of_marg_1(sample[i], p2prob,
			                                        canHaveDuplicates);
		}
	}
	DimT nTg2Dcops = tg2Dcopulas.size();

	double dist = 0.0; // init not really needed - just to avoid a gcc warning..
	double minDist;
	double totDist = 0.0;
	double scProb;

	/// \todo rename bestRows to bestScen
	/// \todo fix the numbering: rows OR scens !!!!

	/*
		The logic in the following is as follows:
		- loop through the ranks of the new margin (sequentially)
		- for each rank, loop through all scenarios and compute the distance
		  (error) of putting the rank to the given scenario,
		  based on the values of already assigned margins
		- store the best scenarios
		- after the loop, chose one of the candidates
	*/
	Copula2D::Cop2DSample::CandList candScens(minNumCandScens, CdfDistEps);
	std::vector<bool> scenUsed(nSc, false);

	std::vector<VectorD> prevRowCdf(nTg2Dcops); // cdf of the previous row
	std::vector<VectorD> colCdfDist(nTg2Dcops); // dist. of using the columns
	for (tg = 0; tg < nTg2Dcops; tg++) {
		prevRowCdf[tg] = ublas::zero_vector<double>(nSc);//Vector(nSc, 0);
		colCdfDist[tg] = VectorD(nSc);
	}

	// compute a safe upper bound for the distance:
	double maxMinDist = 0.0;
	for (tg = 0; tg < nTg2Dcops; tg++) {
		maxMinDist += p2sample2D(oldMargins[tg], marg)->cdfDist(0.0, 1.0);
	}
	maxMinDist *= nSc;

	// loop over ranks of the newly added margin
	for (iR = 0; iR < nSc; iR++) {
		// init minDist & clean the bestRows array
		minDist = maxMinDist;
		candScens.clear();

		for (tg = 0; tg < nTg2Dcops; tg++) {
			// get the dist for the rows
			/// \todo check this!
			p2sample2D(oldMargins[tg], marg)->cdf_dist_of_row(iR, prevRowCdf[tg],
			                                                  colCdfDist[tg]);
		}

		// brute-force approach - this will be SLOW
		for (s = 0; s < nSc; s++) {
			// compute the distance for putting rank iR into scenario s
			if (!scenUsed[s]) {
				dist = 0.0;
				for (tg = 0; tg < nTg2Dcops; tg++) {
					int rankInTgMargin = sample[oldMargins[tg]][s]; // CHECK!
					dist += colCdfDist[tg][rankInTgMargin];
				}
				candScens.insert_cand(s, dist);
			}
		}
		// now we should have a couple..
		assert (candScens.get_num_cand() > 0
		        && "we should have found something...");
		candScens.get_rand_cand(s, minDist);

		TRACE (TrInfo2, "CopulaSample::gen_new_margin(" << marg
			    << "): new link: iR=" << iR << ", s=" << s
			    << "; dist = " << minDist << " (min dist = "
			    << candScens.get_best_cand_value() << ")" << endl);

		scProb = (p2prob == nullptr ? 1.0 / nSc : p2prob[s]);
		totDist += scProb * minDist;
		sample[marg][s] = iR;
		scenUsed[s] = true;

		/// \todo Check!
		for (tg = 0; tg < nTg2Dcops; tg++) {
			i = oldMargins[tg];
			j = sample[i][s]; // the row-rank for the given target copula

			// update prevRowCdf:
			for (DimT jj = j; jj < nSc; jj++) {
				prevRowCdf[tg][jj] += scProb;
			}

			// update p2sample2D
			p2sample2D(i, marg)->add_link(j, iR);

			// update also the [marg][i] copula sample? !!!
		}
	}
	haveSc4Marg[marg] = true;

	// cleaning
	for (i = 0; i < nVar; ++i) {
		for (j = 0; j < nVar; ++j) {
			if (p2sample2D(i, marg) != nullptr) {
				delete p2sample2D(i, marg);
				p2sample2D(i, marg) = nullptr;
				//! temp - should be done more carefully
				// TEST TEMP
				// - delete the whole object, of not shared with other cop. pairs
				// - test that this releases the memory!
				//p2tgInfo(i, marg)->clear_cdf_grid();
				p2tgInfo(i, marg).reset();
			}
		}
	}

	return dist;
}


// the main routine; returns the KS-like-distance
double CopulaSample::gen_sample()
{
	DimT s, marg;
	double totDist = 0.0;

	marg = 0;
	if (!haveSc4Marg[marg]) {
		// initialize the first margin
		for (s = 0; s < nSc; s++) {
			/// \todo Once this works, test it with $nSc - s$ or random numbers.
			/// \todo Make this work with random order, to add the possibility
			///       to shuffle scenarios. An alternative is to shuffle the
			///       output, but it is not so easy as we store the results in
			///       the wrong order (by margin, not by scenario).
			/// \note For now, the values are shuffled at the end,
			///       using CopulaSample::shuffle_results().
			sample[marg][s] = s; // alternative is to use random numbers
		}
		haveSc4Marg[marg] = true;
	}

	MSG (TrInfo, "Starting copula generation.")
	for (marg = 1; marg < nVar; marg++) {
		if (haveSc4Marg[marg]) {
			MSG (TrInfo, " - skipping fixed margin " << marg + 1 << " of " << nVar);
		} else {
			MSG (TrInfo, " - adding margin " << marg + 1 << " of " << nVar);
			totDist += gen_new_margin(marg);
		}
	}
	MSG (TrInfo, "Finished generating the copula.")

	return totDist;
}


// get the ranks as one big matrix
void CopulaSample::get_result_ranks(MatrixI & ranks)
{
	DimT i;
	ranks.resize(nVar, nSc);
	for (i = 0; i < nVar; ++i) {
		ublas::row(ranks, i) = sample[i]; // TEST!
	}
}


void CopulaSample::print_as_txt(string const fName, bool const scaleTo01,
                                int const sortByMarg)
{
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		DimT marg, s;
		if (sortByMarg >= 0) {
			throw std::logic_error(
				"Sorting of output by margins si not yet implemented!");
		}
		//! alt: create a matrix of output values and send it to oFile??
		if (scaleTo01) {
			for (s = 0; s < nSc; s++) {
				for (marg = 0; marg < nVar - 1; marg++) {
					oFile << (static_cast<double>(sample[marg][s]) + 0.5) / nSc << "\t";
				}
				oFile << (static_cast<double>(sample[nVar - 1][s]) + 0.5) / nSc << endl;
			}
		} else {
			// adding 1 to all values, to make ranks go from 1 to nSc
			for (s = 0; s < nSc; s++) {
				for (marg = 0; marg < nVar - 1; marg++) {
					oFile << sample[marg][s] + 1 << "\t";
				}
				oFile << sample[nVar - 1][s] + 1 << endl;
			}
		}
		oFile.close();
	}
}


void CopulaSample::write_gmp_data(string const fName)
{
	DimT i, s;
	std::ofstream oFile;
	oFile.open(fName.c_str(), std::ios::out);
	if (!oFile) {
		cerr << "WARNING: could not open output file " << fName << "!" << endl;
	} else {
		oFile << "# this file was written by the copula-generation code" << endl
		      << endl
		      << "param nVars := " << nVar << ";" << endl
		      << "param nScens := " << nSc << ";" << endl
		      << endl
		      << "param minProbFrac := " << 0.1 << ";" << endl
		      << "param wAvgErr := " << 1.0 << ";" << endl
		      << "param wMaxErr := " << 5.0 << ";" << endl
		      << "param wAvgDev := " << 0.0 << ";" << endl
		      << "param wMaxDev := " << 0.0 << ";" << endl
		      << endl
		      << "param rank (tr)" << endl << "\t:";
		for (i = 0; i < nVar; ++i) {
			oFile << "\t" << i;
		}
		oFile << " :=";
		for (s = 0; s < nSc; ++s) {
			oFile << endl << "\t" << s;
			for (i = 0; i < nVar; ++i) {
				oFile << "\t" << sample[i][s];
			}
		}
		oFile << endl << ";" << endl
		      << endl;

		double tgScProb;
		VectorD uScen(nVar);
		oFile << "param tgCdf :=";
		for (s = 0; s < nSc; ++s) {
			for (i = 0; i < nVar; ++i) {
				uScen(i) = rank2U01(sample[i][s], nSc);
			}
			if (p2copInfo && p2copInfo->has_cdf()) {
				// have the multivar target info
				tgScProb = p2copInfo->cdf(uScen);
			} else {
				// no multivar target -> take the average of all the 2D copulas!
				tgScProb = 0.0;
				DimT nTgCops = 0;
				for (i = 0; i < nVar; ++i) {
					for (DimT j = 0; j < nVar; ++j) {
						if (j != i && p2tgInfo(i, j)) {
							tgScProb += p2tgInfo(i, j)->cdf(uScen(0), uScen(1));
							nTgCops ++;
						}
					}
				}
				tgScProb /= nTgCops;
			}
			oFile << endl << "\t" << s << "\t" << tgScProb;
		}
		oFile << endl << ";" << endl;

		oFile << endl << "end;" << endl;
		oFile.close();
	}
}


void CopulaSample::shuffle_results()
{
	// note: sample is 'std::vector<VectorI>' with dimensions [nVar x nSc]
	DimT i, s;

	// create the vector of shuffled ranks
	VectorI rankOrder(nSc); // 'rankOrder(s) = k' means scen. s gets rank k
	for (s = 0; s < nSc; ++s)
		rankOrder(s) = s;
	std::random_shuffle (rankOrder.begin(), rankOrder.end());

	// we shuffle out-of-place, since it is much easier and less error-prone
	VectorI tmpV;

	for (i = 0; i < nVar; ++i) {
		tmpV = sample[i]; // copy the values
		VectorI & margS = sample[i]; // just a shortcut
		for (s = 0; s < nSc; ++s)
			margS(s) = tmpV(rankOrder(s));
	}
}


void CopulaSample::fix_marg_values(MatrixI const & X, bool const allowDuplicates)
{
	DimT M = X.size1();
	if (M > nVar)
		throw std::length_error("Trying to fix more margins than we have available.");
	if (X.size2() != nSc)
		throw std::length_error("Trying to fix margins with wrong number of scenarios.");

	DimT m;
	for (m = 0; m < M; ++m) {
		if (haveSc4Marg[m])
			throw std::logic_error("Trying to fix an already fixed margin.");
		sample[m] = ublas::row(X, m);
		assert(sample[m].size() == nSc && "dimension check");
		haveSc4Marg[m] = true;
	}

	if (allowDuplicates)
		nmbMargsWithDuplicates = M; // store in the class for later use
}
