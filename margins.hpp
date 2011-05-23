#ifndef MARGINS_HPP
#define MARGINS_HPP

#include "margin-distrib.hpp"


namespace MarginDistrib {

class MarginsInfo {
protected:
	DimT nM; ///< number of margins
	ublas::vector<UnivarMargin::Ptr> p2margins;

public:
	MarginsInfo(DimT const nVars) : nM(nVars), p2margins(nVars) {}
	virtual ~MarginsInfo() {}

	void get_margin_distr(MatrixI const & ranks, MatrixD & values);

	typedef boost::shared_ptr<MarginsInfo> Ptr;
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
	MatrixD & data; // reference, since the data are with the copula class

public:
	SampleMargins(MatrixD & samples,
	              UnivarMargin::SamplePP const postP = UnivarMargin::PP_None);
};

} // namespace MarginDistrib

#endif
