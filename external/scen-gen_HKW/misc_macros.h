#ifndef MISC_MACROS_H
#define MISC_MACROS_H

// ----------------------------------------------------------------
// random-number generation macros

/* Random number generator - double from [0,1]
	 rand() is the only random number generator in C++
	 Note that rand() is from 0 to RAND_MAX inclusive!
	 */
#define frand (((double) rand() + 0.5) / ((double) RAND_MAX + 1.0))

/* Random number generator - integer from 0 to Max-1
	 rand() is the only random number generator in VC++
	 Note that rand() is from 0 to RAND_MAX inclusive!
	 */
#define irand(Max) (rand() % (Max))

/* Random numbers from standard normal distribution
	 using Box-Muller formula
	 */
#define Pi 3.1415926535
#define BM_normal() (double)(sqrt(-2*log(frand))*cos(2*Pi*frand))


// ----------------------------------------------------------------
// gen. mathematical macros

/* min() and max() functions
   - These are defined in MSVC, as long as __STDC__ is undefined
   - With GCC, we can use a special formulation that uses
     GNU C extensions and that avoids double evaluation of f(x)
     in case we cal min(x,f(x)). The `__extension__' keyword
     supresses ISO-C-incompatibility warning with -pedantic.
   - In other cases, we use the standard definitions
*/
#if defined(_MSC_VER) && !defined(__STDC__)
	// defined by compiler -> nothing here
#elif defined(__GNUC__)
	// using GCC extensions
	#define min(X, Y) __extension__       \
	({ typeof (X) __x = (X), __y = (Y);   \
	  (__x < __y) ? __x : __y; })
	#define max(X, Y) __extension__       \
	({ typeof (X) __x = (X), __y = (Y);   \
	  (__x > __y) ? __x : __y; })
#else
	// standard definition
	#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
	#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#endif

// use this for comparing two float numbers
#define EPS 1e-8


// ----------------------------------------------------------------
// generic C/C++ macros

#define allocate(var, type, size) if ((var = (type *) malloc(size * sizeof(type))) == NULL) { printf("\nNOT ENOUGH MEMORY!\n\n");	exit(1); }

#endif
