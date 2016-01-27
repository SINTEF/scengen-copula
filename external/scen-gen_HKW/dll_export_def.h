// ----------------------------------------------------------------
// macros for cross-platform dll-exporting
//
// source: http://gcc.gnu.org/wiki/Visibility

#ifndef HKW_DLL_EXPORT_DEF_H
#define HKW_DLL_EXPORT_DEF_H

#if defined HKW_NO_DLL_DEFS
	// use this if we are just building an executable that does not use DLLs
	#define HKW_DLL_EX
	#define HKW_DLL_IM
#else
	#if defined _WIN32
		#ifdef BUILDING_DLL
			#ifdef __GNUC__
				#define HKW_DLL_EX __attribute__((dllexport))
			#else
				// Note: actually gcc seems to also supports this syntax.
				#define HKW_DLL_EX __declspec(dllexport)
			#endif
		#else
			#ifdef __GNUC__
				#define HKW_DLL_EX __attribute__((dllimport))
			#else
				// Note: actually gcc seems to also supports this syntax.
				#define HKW_DLL_EX __declspec(dllimport)
			#endif
		#endif
		#define HKW_DLL_IM
	#else
		#if __GNUC__ >= 4
			#define HKW_DLL_EX __attribute__ ((visibility("default")))
			#define HKW_DLL_IM __attribute__ ((visibility("hidden")))
		#else
			#define HKW_DLL_EX
			#define HKW_DLL_IM
		#endif
	#endif
#endif

#endif // header guard
