#ifndef BIVARIATE_STUDENT_H
#define BIVARIATE_STUDENT_H


class BivariateCumulativeStudentDistribution
{
private:
	unsigned dof;     ///< degrees of freedom
	double   cor;     ///< correlation
	double const EPS; ///< epsilon for computations

	/// \name parameters used for computation
	///&{
		double unCor;  ///< = 1 - cor^2
		bool dofEven;  ///< is #dof an even number?
	///&}

	/// function x defined on top of page 155
	double f_x(double const m, double const h, double const k) const;

	double P_n(double const h, double const k) const;

protected:

public:
	/** Default constructor */
	BivariateCumulativeStudentDistribution(unsigned const degF, double const rho);
	/** Default destructor */
	~BivariateCumulativeStudentDistribution();

	double cdf(double const x, double const y) const { return P_n(x, y); }

	double operator() (double const a, double const b) const
	{ return P_n(a, b); }
};

#endif // BIVARIATE_STUDENT_H
