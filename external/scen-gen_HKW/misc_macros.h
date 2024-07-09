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

/* min() and max() macros
   With the 'standard definitions' (see below),
   min(x,f(x)) will cause f(x) to be computed twice.
   This can be avoided using GNU C extensions:
*/
#ifdef __GNUC__
	#define min(X, Y)            \
	({                           \
		__typeof__ (X) _x = (X); \
		__typeof__ (Y) _y = (Y); \
		_x < _y ? _x : _y;       \
	})
	#define max(X, Y)            \
	({                           \
		__typeof__ (X) _x = (X); \
		__typeof__ (Y) _y = (Y); \
		_x > _y ? _x : _y;       \
	})
#else
	// Standard definitions:
	#ifndef min
		#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
		#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
	#endif
#endif

// use this for comparing two float numbers
#define EPS 1e-8


// ----------------------------------------------------------------
// generic C/C++ macros

#define allocate(var, type, size) if ((var = (type *) malloc(size * sizeof(type))) == NULL) { printf("\nNOT ENOUGH MEMORY!\n\n");	exit(1); }

#endif // header guard
