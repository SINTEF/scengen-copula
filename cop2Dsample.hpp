#ifndef COP_2D_SAMPLE_HPP
#define COP_2D_SAMPLE_HPP

#include <string>
//#include <cassert>

#include "common.hpp"
#include "cop2Dinfo.hpp"

using std::string;


// forward declarations
// needed to make CopulaSample friend of Cop2DSample
namespace CopulaScen {
	class CopulaSample;
}


namespace Copula2D {


/// sample from a bivariate copula in terms of ranks
/**
	Internally, the class only works on connections between ranks.
	To be able to connect this to scenarios, we keep extra vectors to convert
	vector indices to scenario numbers.
**/
class Cop2DSample {
	friend class CopulaScen::CopulaSample;
private:
	/// \name main parameters defining the copula
	//@{
		int N;        ///< number of sample points in the copula
		IVector i2jC; ///< vector of i->j connections
		IVector j2iC; ///< vector of j->i connections

		void update_j2iC_from_i2jC(); ///< synchronize j2iC with i2jC
		void update_i2jC_from_j2iC(); ///< synchronize i2jC with j2iC

		/// check consistency of the i->j and j->i links
		bool have_valid_cop();
	//@}

	/// \name matrix with the target values of the cdf - speeds things up
	///@{
		Matrix tgCdfOfR;      ///< pre-computed matrix of target cdf values
		void fill_tgCdfOfR(); ///< get the matrix of target cdf values
	///@}

	/// \name connection to scenarios of the main alg.
	//@{
		/// pointer to scenario probabilities.
		/// can be NULL, implying equiprobable scenarios/samples
		double const *p2prob;

		/// scenarios of ranks of the 1st margin.
		/// use \c p2prob[scenOfMarg1R[i]] to get probability of rank \c i
		IVector scenOfMarg1R;

		//Vector cumProb;     ///< cummulative probabilities

		/// position of the cdf evaluation point in the discretization intervals.
		/// reasonable values are 0.5 (cond. mean) or 1.0 (standard emp. cdf)
		double evalPtPos;

		/// discretization points of the copula marginals.
		/// these are sorted as 'i', use scenOfRow[] to get scenario values
		Vector copEvalPts;

		string sampleId; ///< string to identify the sample (for reporting)
	//@}

	Cop2DInfo const *p2tgInfo; ///< pointer to the target specification (can be 0)

	/// \name misc. methods
	//@{
		/// compute probability (nmb. of pts) in a given box
		double prob_in_box(int const i0, int const j0,
										   int const i1, int const j1) const;

		/// convert ranks to the unif[0,1) values
		inline double rank2U01(int const r) const {
			return static_cast<double>(r) / N;
		}

		/// convert unif[0,1) values to ranks
		inline int u01ToRank(double const z) const { return lround(N * z); }

		/// convert ranks to the copula margins, i.e. val in (0,1)
		//inline double rank2cop(int const i) const { return (i + 0.5) / N; }

		/// convert copula margins, i.e. val in (0,1), to ranks
		//inline int cop2rank(double const z) const { return lround(N * z - 0.5); }

		/// shortcut to the target rank-cdf, with margins given as ranks
		inline double tgRCdfOfR(int const i, int const j) const {
			return N * tgCdfOfR[i][j];
		}

		void gen_random();
	//@}

	/// distance used when comparing cdf to the target (point-wise).
	/// we can use for ex. square to get a bigger penalty on big deviations
	/// u and v should be from (0, 1)
	inline double cdfDist(double const u, double const v) const {
		assert (u >= 0.0 && u <= 1.0 + DblEps && v >= 0.0 && v <= 1.0 + DblEps
						&& "bound check");
		return fabs(u - v);
	}

	/// distance used when comparing rank-cdf to the target (point-wise).
	/// uR and vR should be from (0, N)
	inline double cdfDistOfR(int const uR, int const vR) const {
		return cdfDist(rank2U01(uR), rank2U01(vR));
	}

protected:

public:
	Cop2DSample(int const nSamples, Cop2DInfo const *const p2TgCop,
							string const id = "");

	~Cop2DSample() {};

	/// set the scenarios for the 'i' variable
	/**
		\param[in] scenOfColR vector of scenarios for the ranks of column \c i
		           - \code scenOfColR[i] = s \endcode means that scenario \c s
		             includes (consists of) rank \c i in the first variable
		             of the bivariate sample (i.e. column \c i).
		\param[in] p2scProb pointer to scenario probabilities;
		           NULL means equiprobable scenarios
	**/
	//void set_scen_of_i(IVector const &scenOfColR, double const *p2scProb = NULL);

	/// set the scenarios for the the first margin
	/**
		\param[in] margScen vector of scenarios for the margin
		           - ordered by scenarios, i.e. \code margScen[s] = i \endcode
		           means that scenario \c s includes (consists of) rank \c i in
		           the first variable of the bivariate sample.
		\param[in] p2scProb pointer to scenario probabilities;
		           NULL means equiprobable scenarios
	**/
	void set_scen_of_marg_1(IVector const &margScen,
													double const *p2scProb = NULL);

	/// the main heuristics
	double gen_heur();

	/// cdf -> returns values between 0 and 1
	double cdfOfR(int const i, int const j) const;

	/// rank cdf -> returns integral values between 0 and N-1
	/// makes sense only in the equiprobable case!
	int rCdfOfR(int const i, int const j) const {
		return u01ToRank(cdfOfR(i, j));
	}

	void print_as_txt(string const fName, bool const sortByScen = false);

	/// Compute the cdf-distance of for a whole column, given Cdf of prev. col.
	/**
		\param[in]  i index of the column, i.e. the rank \c i
		\param[in]  prevColCdf vector of \f$ F(i-1,\cdot) \f$ values
		\param[out] rowCdfDist cdf-distance for putting the link to different rows
		\param[out] rowFree if given, will be filled by indicators whether the
		            rows are available or not; in this case, the rowCdfDist for
		            the rows with 'false' need not be defined.
	**/
	void cdf_dist_of_col(int const i, Vector const &prevColCdf,
											 Vector &rowCdfDist, bool rowFree[] = NULL);

	/// Compute the cdf-distance of for a whole row, given Cdf of prev. row.
	/**
		\param[in]  j index of the row, i.e. the rank \c j
		\param[in]  prevRowCdf vector of \f$ F(\cdot,j-1) \f$ values
		\param[out] colCdfDist cdf-distance for putting the link to different cols
		\param[out] colFree if given, will be filled by indicators whether the
		            columns are available or not; in this case, the colCdfDist for
		            the columns with 'false' need not be defined.
	**/
	void cdf_dist_of_row(int const j, Vector const &prevRowCdf,
											 Vector &colCdfDist, bool colFree[] = NULL);

	/// \todo do this; return false if this breaks another link
	bool add_link(int const colR, int const rowR);
};


} // namespace

#endif
