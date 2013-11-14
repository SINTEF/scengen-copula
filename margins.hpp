#ifndef MARGINS_HPP
#define MARGINS_HPP

#include "margin-distrib.hpp"


namespace MarginDistrib {

class MarginsInfo {
protected:
	DimT nM; ///< number of margins
	/// pointers to the univariate margins (they provide CDF etc)
	ublas::vector<UnivarMargin::Ptr> p2margins;

public:
	MarginsInfo(DimT const nVars) : nM(nVars), p2margins(nVars) {}
	MarginsInfo(std::string fName);
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
	void read_from_file(std::string fName);
public:
	MixedMargins(DimT const nVars) : MarginsInfo(nVars) {}
	MixedMargins(std::string fName);

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
