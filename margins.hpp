#ifndef MARGINS_HPP
#define MARGINS_HPP

#include "margin-distrib.hpp"

#include <map>


namespace MarginDistrib {

/// \name objects for the margin-type name map, used in the main code
///@{
	/// enum for the known multivariate margin specification types
	enum class MargTypeID {normal, sample, fixed, mixed, unknown};

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
	DistribNameMapT distrNameMap;

public:
	MixedMargins(DimT const nVars) : MarginsInfo(nVars) {}

	/// constructor based on a file
	/**
		\param[in] tgFName  name of the input file
		\param[in]   nVars  number of variables - optional, used for checking
	**/
	MixedMargins(std::string const & tgFName, DimT const nVars = 0);

	/// attach margin
	void attach_margin(UnivarMargin::Ptr & p2marg, DimT const index);
};


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


class SampleMargins : public MarginsInfo {
private:
	void init_margins(UnivarMargin::SamplePP const postP);
	MatrixD const & data; // reference, since the data are with the copula class

public:
	SampleMargins(MatrixD const & samples,
	              UnivarMargin::SamplePP const postP = UnivarMargin::PP_None);
};

} // namespace MarginDistrib

#endif
