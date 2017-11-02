#ifndef HVL_H
#define HVL_H

#include <hvl_config.h>

#ifdef __cplusplus
	#include <cstddef>
#else
	#include <stddef.h>
#endif // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifdef LIBHVL_PLATFORM_WIN32
	#define LIBHVL_DLLCALL __cdecl

	#ifdef LIBHVL_DLLBUILD
		#if defined LIBHVL_COMPILER_MINGW || defined LIBHVL_COMPILER_MSYS
			#define LIBHVL_DLLEXPORT __attribute__((dllexport))
		#elif defined LIBHVL_COMPILER_MSVC
			#define LIBHVL_DLLEXPORT __declspec(dllexport)
		#else
			#error "Unsupported or misdetected compiler"
		#endif // LIBHVL_COMPILER_
	#else
		#if defined LIBHVL_COMPILER_MINGW || defined LIBHVL_COMPILER_MSYS
			#define LIBHVL_DLLEXPORT __attribute__((dllexport))
		#elif defined LIBHVL_COMPILER_MSVC
			#define LIBHVL_DLLEXPORT __declspec(dllexport)
		#else
			#error "Unsupported or misdetected compiler"
		#endif // LIBHVL_COMPILER_
	#endif // LIBHVL_DLLBUILD
#elif defined LIBHVL_PLATFORM_UNIX
	#ifndef __x86_64__
		#define LIBHVL_DLLCALL __attribute__((__cdecl__))
	#else
		#define LIBHVL_DLLCALL
	#endif // __x86_64__

	#define LIBHVL_DLLEXPORT __attribute__((visibility("default")))
	#define LIBHVL_THREADING_PTHREAD
#endif // LIBHVL_PLATFORM_


struct HVL_Context;
struct HVL_Prepared;

typedef struct HVL_Context HVL_Context;
typedef struct HVL_Prepared HVL_Prepared;

typedef struct HVL_Pair {
	double x;
	double y;
} HVL_Pair;
typedef struct HVL_Range {
	HVL_Pair *p;
	size_t count;
} HVL_Range;

typedef enum HVL_RetCode {
	HVL_OK,
	HVL_E_NO_MEMORY,
	HVL_E_INVALID_ARG
} HVL_RetCode;

LIBHVL_DLLEXPORT void LIBHVL_DLLCALL HVL_free_context(HVL_Context *ctx);
LIBHVL_DLLEXPORT void LIBHVL_DLLCALL HVL_free_prepared(HVL_Prepared *prep);
LIBHVL_DLLEXPORT void LIBHVL_DLLCALL HVL_free_range(HVL_Range *r);
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL HVL_make_context(HVL_Context **ctx, const int prec);
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL HVL_prepare(HVL_Prepared **prep, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_t_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_dx_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da0_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da1_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da2_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da3_prepared(const HVL_Prepared *prep);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_t(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_dx(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da0(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da1(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da2(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT double LIBHVL_DLLCALL HVL_da3(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_t_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_dx_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_da0_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_da1_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_da2_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL HVL_da3_range(const HVL_Context *ctx, const double from, const double to, const double step, const double a0, const double a1, const double a2, const double a3);

#ifdef __cplusplus
}

#include <type_traits>
static_assert(std::is_pod<HVL_Pair>::value, "HVL_Pair is not a POD type");
static_assert(std::is_pod<HVL_Range>::value, "HVL_Range is not a POD type");

#endif // __cplusplus

#endif // HVL_H
