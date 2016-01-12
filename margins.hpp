#ifndef MARGINS_HPP
#define MARGINS_HPP

#include "margin-distrib.hpp"

#include <map>


namespace MarginDistrib {

/// \name objects for the margin-type name map, used in the main code
///@{
	/// enum for the known multivariate margin specification types
	enum class MargTypeID {normal, sample, moments, fixed, mixed, unknown};

	/// type for the copula map
	typedef std::map<std::string, MargTypeID> MargNameMapT;

	/// this fills the copula map with the known copula types
	void make_marg_name_map(MargNameMapT & mMap);
///@}


class MarginsInfo {
protected:
	DimT nM; ///< number of margins
	/// pointers to the univariate margins (they provide CDF etc)
	ublas::vector<UnivarMargin::Ptr> p2margins;

public:
	MarginsInfo(DimT const nVars) : nM(nVars), p2margins(nVars) {}
	//MarginsInfo(std::string fName);
	virtual ~MarginsInfo() {}

	/// convert ranks to distribution values
	/**
		\param[in]   ranks  matrix of ranks [nMargins * nScens]
		\param[out] values  matrix of values [nMargins * nScens]
	**/
	void get_margin_distr(MatrixI const & ranks, MatrixD & values);

	/// get dimension
	DimT dim() const { return nM; }

	typedef boost::shared_ptr<MarginsInfo> Ptr;
};


class MixedMargins : public MarginsInfo {
private:
	//DistribNameMapT distrNameMap;

public:
	MixedMargins(DimT const nVars) : MarginsInfo(nVars) {}

	/// constructor based on a file
	/**
		\param[in] tgFName  name of the input file
		\param[in]  nScens  number of scenarios - optional, but needed for moments
		\param[in]   nVars  number of variables - optional, used for checking
	**/
	MixedMargins(std::string const & tgFName, DimT const nScens = 0,
	             DimT const nVars = 0);

	/// attach margin
	void attach_margin(UnivarMargin::Ptr & p2marg, DimT const index);
};


#ifdef HAS_QUANTLIB
class NormalMargins : public MarginsInfo {
private:
	VectorD means;  ///< vector of means
	VectorD stDevs; ///< vector of standard deviations

	void read_from_file(std::string fName);

	void init_margins(UnivarMargin::SamplePP const postP);

public:
	NormalMargins(VectorD const & EVs, VectorD const & SDs,
	              UnivarMargin::SamplePP const postP = UnivarMargin::PP_None);

	NormalMargins(std::string fName,
	              UnivarMargin::SamplePP const postP = UnivarMargin::PP_None);
};
#endif // HAS_QUANTLIB

class SampleMargins : public MarginsInfo {
private:
	void init_margins(UnivarMargin::SamplePP const postP);
	MatrixD const & data; // reference, since the data are with the copula class

public:
	/// constructor with the target samples (matrix is [nVar * nPts])
	SampleMargins(MatrixD const & samples,
	              UnivarMargin::SamplePP const postP = UnivarMargin::PP_None);
};


#ifdef HAS_HKW
/// class for all margins specified by their moments
/**
	\warning needs the moment-matching library by Michal Kaut

	\todo make it a class derived from fixed margins?

	\note this is just a wrapper, the actual work is done by the univariate class
**/
class MarginsByMoments : public MarginsInfo {
private:
	std::vector<VectorD> tgMoments; ///< target moments
	int formOfMoments;              ///< format of moments, values as in HKW
	DimT nSc;                       ///< number of scenarios

	/// read the target moments from a file
	/**
		\param[in]  fName  name of the file with the moments
		\param[in] matFmt  does the input file use the (old) matrix-style format?
	**/
	void read_from_file(std::string fName, bool const matFmt = false);

	/// initialize the univariate classes (generate scenarios)
	void init_margins();

public:
	/// constructor with file name and other parameters
	/**
		\param[in]  fName  name of the file with the moments
		\param[in]  nScen  number of scenarios to generate
		\param[in]    FoM  format of moments in the input files
		\param[in] matFmt  does the input file use the (old) matrix-style format?
	**/
	MarginsByMoments(std::string fName, DimT const nScen, int const FoM = 0,
	                 bool const matFmt = false);
};
#endif // HAS_HKW

} // namespace MarginDistrib

#endif
