#include "bivariate_student.h"

#include <cmath>
#include <cassert>

double const PI = 3.1415926535897932384626433832795;
double const TWOPI = 2 * PI;


// --------------------------------------------------------------------------
// functions used to compute the bivariate cdf, local to this file

/// signum (x)
/**
	taken from http://stackoverflow.com/a/4609795
**/
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

/// arctan of a fraction
/**
	The paper uses scales this to [0,2 pi], unlike the atan2 function in C++
	that gives results in [-pi,pi].
**/
double atanFrac(double const x, double const y) {
	double res = std::atan2(x, y);
	if (res < 0)
		res += 2 * PI;
	assert (
		(  (x >= 0 && y >= 0 && res >= 0 && res <= PI / 2)
		|| (x >= 0 && y <= 0 && res >= PI / 2 && res <= PI)
		|| (x <= 0 && y <= 0 && res >= PI && res <= 1.5 * PI)
		|| (x <= 0 && y >= 0 && res >= 1.5 * PI && res <= TWOPI)
		) && "range check" );
	return res;
}


// --------------------------------------------------------------------------
// class BivariateCumulativeStudentDistribution

// function x(m,h,k) defined on top of page 155
double BivariateCumulativeStudentDistribution::f_x(
	double const m, double const h, double const k) const
{
	double sub = pow(h - cor * k, 2);
	double denom = sub + unCor * (m + pow(k, 2));
	if (denom < EPS)
		return 0; // limit case for cor = +/-1.0
	double res = sub / (sub + unCor * (m + pow(k, 2)));
	assert (res >= 0 && res <= 1 && "range check");
	return res;
}


// this calculates the cdf
double BivariateCumulativeStudentDistribution::P_n(
	double const h, double const k) const
{
	unsigned n = dof; // use the notation from the paper
	double rho = cor; // use the notation from the paper

	unsigned j;
	double res;

	double div = 4 * sqrt(n * PI);
	double xHK = f_x(n, h, k);
	double xKH = f_x(n, k, h);
	double divH = 1 + pow(h, 2) / n;
	double divK = 1 + pow(k, 2) / n;
	double sgnHK = sgn(h - rho * k);
	double sgnKH = sgn(k - rho * h);

	if (dofEven) { // equation (10)
		// first line of (10)
		res = atanFrac(sqrt(unCor), -rho) / TWOPI;

		// second line of (10)
		double dgM = 2 * (1 - xHK);  // multiplier for dgj
		double gjM = sgnHK * 2 / PI; // multiplier for g_j
		// initializations for j = 1:
		double f_j = sqrt(PI / divK);
		double g_j = 1 + gjM * atanFrac(sqrt(xHK), sqrt(1 - xHK));
		double sum = f_j * g_j;
		if (n >= 4) {
			// different formulas for j = 2:
			f_j *= 0.5 / divK; // (2 - 1.5) / (double) (2 - 1) / divK;
			double dgj = gjM * sqrt(xHK * (1 - xHK));
			g_j += dgj;
			sum += f_j * g_j;
			// and then the loop for the rest of the j's:
			for (j = 3; j <= n / 2; ++j) {
				f_j *= (j - 1.5) / (double) (j - 1) / divK;
				dgj *= (double) (j - 2) / (2 * j - 3) * dgM;
				g_j += dgj;
				sum += f_j * g_j;
			}
		}
		res += k / div * sum;

		// third line of (10)
		dgM = 2 * (1 - xKH);
		gjM = sgnKH * 2 / PI;
		// initializations for j = 1:
		f_j = sqrt(PI / divH);
		g_j = 1 + gjM * atanFrac(sqrt(xKH), sqrt(1 - xKH));
		sum = f_j * g_j;
		if (n >= 4) {
			// different formulas for j = 2:
			f_j *= 0.5 / divH; // (2 - 1.5) / (double) (2 - 1) / divK;
			double dgj = gjM * sqrt(xKH * (1 - xKH));
			g_j += dgj;
			sum += f_j * g_j;
			// and then the loop for the rest of the j's:
			for (j = 3; j <= n / 2; ++j) {
				f_j *= (j - 1.5) / (double) (j - 1) / divH;
				dgj *= (double) (j - 2) / (2 * j - 3) * dgM;
				g_j += dgj;
				sum += f_j * g_j;
			}
		}
		res += h / div * sum;

	} else { // equation (11)
		// first line of (11)
		double hk = h * k;
		double hkcn = hk + rho * n;
		double sqrtExpr = sqrt(pow(h, 2) - 2 * rho * hk + pow(k, 2) + n * unCor);
		res = atanFrac(sqrt(n) * (-(h + k) * hkcn - (hk - n) * sqrtExpr),
		               (hk - n) * hkcn - n * (h + k) * sqrtExpr ) / TWOPI;

		if (n > 1) {
			// second line of (11)
			double mult = (1 - xHK) / 2;
			// initializations for j = 1:
			double f_j = 2 / sqrt(PI) / divK;
			double dgj = sgnHK * sqrt(xHK);
			double g_j = 1 + dgj;
			double sum = f_j * g_j;
			// and then the loop for the rest of the j's:
			for (j = 2; j <= (n - 1) / 2; ++j) {
				f_j *= (double) (j - 1) / (j - 0.5) / divK;
				dgj *= (double) (2 * j - 3) / (j - 1) * mult;
				g_j += dgj;
				sum += f_j * g_j;
			}
			res += k / div * sum;

			// third line of (11)
			mult = (1 - xKH) / 2;
			// initializations for j = 1:
			f_j = 2 / sqrt(PI) / divH;
			dgj = sgnKH * sqrt(xKH);
			g_j = 1 + dgj;
			sum = f_j * g_j;
			// and then the loop for the rest of the j's:
			for (j = 2; j <= (n - 1) / 2; ++j) {
				f_j *= (double) (j - 1) / (j - 0.5) / divH;
				dgj *= (double) (2 * j - 3) / (j - 1) * mult;
				g_j += dgj;
				sum += f_j * g_j;
			}
			res += h / div * sum;
		}
	}
	return res;
}


BivariateCumulativeStudentDistribution::BivariateCumulativeStudentDistribution(
	unsigned const degF, double const rho)
: dof(degF), cor(rho), EPS(1e-8),
  unCor(1 - pow(rho, 2)), dofEven(degF % 2 == 0)
{
	//ctor
}

BivariateCumulativeStudentDistribution::~BivariateCumulativeStudentDistribution()
{
	//dtor
}
