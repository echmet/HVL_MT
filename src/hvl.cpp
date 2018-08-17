#include <hvl.h>
#include <new>
#include <cstdlib>

/* Use appropriate headers */
#ifdef LIBHVL_USE_MPFR
	#include "mpreal.h"
#else
	#define _USE_MATH_DEFINES
	#include <cmath>
#endif // HVL_USE_MPFR

#ifdef LIBHVL_PLATFORM_WIN32
	#include <windows.h>
#elif defined LIBHVL_PLATFORM_UNIX
	#include <sys/sysinfo.h>
#endif // LIBHVL_PLATFORM_

#ifdef LIBHVL_THREADING_PTHREAD
	#include <pthread.h>
#endif // LIBHVL_THREADING_PTHREAD

/* Define float variable to appropriate type */
#ifdef LIBHVL_USE_MPFR
	typedef mpfr::mpreal hvl_float;
#else
	typedef double hvl_float;
#endif // HVL_USE_MPFR

/* Provide convertors between native hvl_float type and machine double */
#ifdef LIBHVL_USE_MPFR
	#define ToDouble(x) (x).toDouble()
	#define ToHVLFloat(x) mpfr::mpreal(x)
#else
	#define ToDouble(x) static_cast<double>(x)
	#define ToHVLFloat(x) x
#endif // HVL_USE_MPFR

/* Internal library functions */

/** Internal math */

inline hvl_float Exp(const hvl_float &x)
{
#ifdef LIBHVL_USE_MPFR
	return mpfr::exp(x);
#else
	return std::exp(x);
#endif
}

inline hvl_float Erf(const hvl_float &x)
{
#ifdef LIBHVL_USE_MPFR
	return mpfr::erf(x);
#else
	return std::erf(x);
#endif
}

inline hvl_float Pow(const hvl_float &base, const hvl_float &exp)
{
#ifdef LIBHVL_USE_MPFR
	return mpfr::pow(base, exp);
#else
	return std::pow(base, exp);
#endif
}

inline hvl_float Sqrt(const hvl_float &x)
{
#ifdef LIBHVL_USE_MPFR
	return mpfr::sqrt(x);
#else
	return std::sqrt(x);
#endif
}

#ifdef LIBHVL_USE_MPFR

#include <csignal>
#include <csetjmp>
#include <mutex>
#include <time.h>

#ifdef _WIN32
#define THREAD_SLEEP(n) Sleep(n)
#else
#include <unistd.h>
#define THREAD_SLEEP(n) sleep(n)
#endif // _WIN32


#ifdef LIBHVL_THREADING_PTHREAD
	typedef pthread_t thread_id_t;
	#define CURRENT_THREAD_ID pthread_self
#elif defined LIBHVL_THREADING_WIN32
	typedef DWORD thread_id_t;
	#define CURRENT_THREAD_ID GetCurrentThreadId
#endif // LIBHVL_THREADING_

#ifdef LIBHVL_CATCH_MPFR_ASSERTS

static jmp_buf env;
static thread_id_t mpfr_failure_parent_thread_id;
static std::mutex mpfr_assert_handler_mtx;
#ifndef _WIN32
typedef void (*abrthandler_t)(int, siginfo_t *, void *);

struct sigaction mpfr_crashAction;
struct sigaction mpfr_origAction;
#else
typedef _invalid_parameter_handler abrthandler_t;

static abrthandler_t mpfr_prev_handler;
#endif // _WIN32
static volatile bool mpfr_assert_triggered;

static
#ifdef _WIN32
void __cdecl mpfr_assertion_failed_handler(const wchar_t *, const wchar_t *, const wchar_t *, unsigned int, uintptr_t)
#else
void mpfr_assertion_failed_handler(int, siginfo_t *, void *)
#endif // _WIN32
{
	mpfr_assert_handler_mtx.lock();

	if (mpfr_failure_parent_thread_id != CURRENT_THREAD_ID()) {
		mpfr_assert_triggered = true;
		mpfr_assert_handler_mtx.unlock();
		for (;;)
			THREAD_SLEEP(1);
	}

	if (mpfr_assert_triggered) {
		mpfr_assert_handler_mtx.unlock();
		THREAD_SLEEP(5);
		return;
	}


	mpfr_assert_triggered = true;
	mpfr_assert_handler_mtx.unlock();

	longjmp(env, 1);
}


#ifdef _WIN32
	static
	void install_abort_handler(abrthandler_t handler)
	{
		mpfr_prev_handler = _set_invalid_parameter_handler(handler);
	}

	static
	void restore_abort_handler()
	{
		_set_invalid_parameter_handler(mpfr_prev_handler);
	}
#else
	static
	void install_abort_handler(abrthandler_t handler)
	{
		memset(&mpfr_crashAction, 0, sizeof(struct sigaction));
		memset(&mpfr_origAction, 0, sizeof(struct sigaction));
		sigemptyset(&mpfr_crashAction.sa_mask);
		sigaddset(&mpfr_crashAction.sa_mask, SIGABRT);

		mpfr_crashAction.sa_flags = SA_SIGINFO;
		mpfr_crashAction.sa_sigaction = handler;
		sigaction(SIGABRT, &mpfr_crashAction, &mpfr_origAction);
	}

	static
	void restore_abort_handler()
	{
		sigaction(SIGABRT, &mpfr_origAction, nullptr);
	}
#endif // _WIN32

#define MPFR_TRY_BEG \
	if (!setjmp(env)) { \
		mpfr_failure_parent_thread_id = CURRENT_THREAD_ID(); \
		install_abort_handler(mpfr_assertion_failed_handler); \
		mpfr_assert_triggered = false;

#define MPFR_TRY_END \
		restore_abort_handler(); \
	}

#define MPFR_CATCH_BEG \
	else {
#define MPFR_CATCH_END \
	}

#else

#define MPFR_TRY_BEG
#define MPFR_TRY_END
#define MPFR_CATCH_BEG
#define MPFR_CATCH_END

static bool mpfr_assert_triggered{false};

#endif  // LIBHVL_CATCH_MPFR_ASSERTS

#endif // LIBHVL_USE_MPFR

/** Internal data types */
struct HVL_Context {
	explicit HVL_Context(const int prec_bits) :
	#ifdef LIBHVL_USE_MPFR
		HVL_SQRT2(Sqrt(mpfr::mpreal("2.0"))),
		HVL_PI(mpfr::const_pi()),
		HVL_SQRTPI(Sqrt(mpfr::const_pi())),
		HVL_ZERO(mpfr::mpreal("0.0")),
		HVL_HALF(mpfr::mpreal("0.5")),
		HVL_ONE(mpfr::mpreal("1.0")),
		HVL_TWO(mpfr::mpreal("2.0")),
		HVL_SQRT2REC(mpfr::mpreal("1.0") / Sqrt(mpfr::mpreal("2.0"))),
	#else
		HVL_SQRT2(Sqrt(2.0)),
		HVL_PI(M_PI),
		HVL_SQRTPI(Sqrt(M_PI)),
		HVL_ZERO(0.0),
		HVL_HALF(0.5),
		HVL_ONE(1.0),
		HVL_TWO(2.0),
		HVL_SQRT2REC(1.0 / Sqrt(2.0)),
	#endif // HVL_USE_MPFR
	prec_bits(prec_bits)
	{}

	explicit HVL_Context(const HVL_Context &other) :
		HVL_SQRT2(other.HVL_SQRT2),
		HVL_PI(other.HVL_PI),
		HVL_SQRTPI(other.HVL_SQRTPI),
		HVL_ZERO(other.HVL_ZERO),
		HVL_HALF(other.HVL_HALF),
		HVL_ONE(other.HVL_ONE),
		HVL_TWO(other.HVL_TWO),
		HVL_SQRT2REC(other.HVL_SQRT2REC),
		prec_bits(other.prec_bits)
	{}

	const hvl_float HVL_SQRT2;
	const hvl_float HVL_PI;
	const hvl_float HVL_SQRTPI;
	const hvl_float HVL_ZERO;
	const hvl_float HVL_HALF;
	const hvl_float HVL_ONE;
	const hvl_float HVL_TWO;
	const hvl_float HVL_SQRT2REC;
	const long int prec_bits;
};

struct HVL_Prepared {
	HVL_Prepared(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3) :
		ctx(HVL_Context(*ctx)),
		x(x),
		a0(a0),
		a1(a1),
		a2(a2),
		a3(a3),
		sqrtz(0), nmr(0), gauss_base(0), b1pb2(0), reca3(0), xmu(0) /* These will be set to correct values externally later */
	{}

	const HVL_Context ctx;

	/* HVL function parameters */
	const hvl_float x;
	const hvl_float a0;
	const hvl_float a1;
	const hvl_float a2;
	const hvl_float a3;

	/* Precomputed values used to calculate HVL and all of its derivatives */
	const hvl_float sqrtz;		/* Used directly by da3 only */
	const hvl_float nmr;		/* Used by dx, da1 and da2 */
	const hvl_float gauss_base;
	const hvl_float b1pb2;
	const hvl_float reca3;
	const hvl_float xmu;
};

/** Intermediate results calculators */
inline hvl_float CalcGaussBase(const HVL_Context *ctx, const hvl_float &q, const hvl_float &expSqrtz2)
{
	return expSqrtz2 * q / ctx->HVL_SQRTPI;
}

inline hvl_float CalcB1pB2(const HVL_Context *ctx, const hvl_float &a3, const hvl_float &sqrtz)
{
	const hvl_float b1 = ctx->HVL_ONE / (Exp(a3) - ctx->HVL_ONE);
	const hvl_float b2 = ctx->HVL_HALF * (ctx->HVL_ONE + Erf(sqrtz));

	return b1 + b2;
}

inline hvl_float CalcNmr(const HVL_Context *ctx, const hvl_float &a0, const hvl_float &q, const hvl_float &expSqrtz2)
{
	return a0 * Pow(q, 2) * Pow(expSqrtz2, 2) / ctx->HVL_PI;
}

inline hvl_float CalcReca3(const HVL_Context *ctx, const hvl_float &a3)
{
	return ctx->HVL_ONE / a3;
}

inline hvl_float CalcXmu(const hvl_float &x, const hvl_float &a1)
{
	return x - a1;
}

inline hvl_float CalcExpSqrtz2(const hvl_float &sqrtz)
{
	return Exp(-sqrtz*sqrtz);
}

inline hvl_float CalcQ(const HVL_Context *ctx, const hvl_float &a2)
{
	return ctx->HVL_SQRT2REC / a2;
}

inline hvl_float CalcSqrtz(const hvl_float &xmu, const hvl_float &q)
{
	return xmu * q;
}

inline void CheckPrecision(const HVL_Context *ctx)
{
#ifdef LIBHVL_USE_MPFR
	const long int def_prec = mpfr::mpreal::get_default_prec();

	if (ctx->prec_bits != def_prec)
		mpfr::mpreal::set_default_prec(ctx->prec_bits);
#else
	(void)ctx; // No-op
#endif // LIBHVL_USE_MPFR
}

/** Internal HVL calculators which use the intermediate results instead of the basic HVL parameters */
inline double HVL_t_internal(const HVL_Context *ctx, const hvl_float &a0, const hvl_float &a3, const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3)
{
	const hvl_float gauss = a0 * gauss_base;

	if (a3 == ctx->HVL_ZERO)
		return ToDouble(gauss);

	return ToDouble(reca3 * gauss / b1pb2);
}

inline double HVL_dx_internal(const HVL_Context *ctx,
			      const hvl_float &a0, const hvl_float &a2, const hvl_float &a3,
			      const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3, const hvl_float &nmr, const hvl_float &xmu)
{
	const hvl_float gauss_der = gauss_base * -a0 * xmu / (a2 * a2);

	if (a3 == ctx->HVL_ZERO)
		return ToDouble(gauss_der);

	const hvl_float result = reca3 * ((gauss_der / b1pb2) -nmr / (b1pb2 * b1pb2));

	return ToDouble(result);

}

inline double HVL_da0_internal(const HVL_Context *ctx,
			       const hvl_float &a3,
			       const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3)
{
	if (a3 == ctx->HVL_ZERO)
		return ToDouble(gauss_base);

	return ToDouble(reca3 * gauss_base / b1pb2);
}

inline double HVL_da1_internal(const HVL_Context *ctx,
			       const hvl_float &a0, const hvl_float &a2, const hvl_float &a3,
			       const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3, const hvl_float &nmr, const hvl_float &xmu)
{
	const hvl_float gauss_der = a0 * xmu / (a2 * a2) * gauss_base;

	if (a3 == ctx->HVL_ZERO)
		return ToDouble(gauss_der);

	const hvl_float result = reca3 * ((gauss_der / b1pb2) + nmr / (b1pb2 * b1pb2));

	return ToDouble(result);
}

inline double HVL_da2_internal(const HVL_Context *ctx,
			       const hvl_float &a0, const hvl_float &a2, const hvl_float &a3,
			       const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3, const hvl_float &nmr, const hvl_float &xmu)
{
	const hvl_float gauss_der = a0 * (xmu * xmu * gauss_base / Pow(a2, 3) - gauss_base / a2);

	if (a3 == ctx->HVL_ZERO)
		return ToDouble(gauss_der);

	const hvl_float rnmr = nmr * xmu / a2;
	const hvl_float result = reca3 * ((gauss_der / b1pb2) + rnmr / (b1pb2 * b1pb2));

	return ToDouble(result);
}

inline double HVL_da3_internal(const HVL_Context *ctx,
			       const hvl_float &a0, const hvl_float &a3,
			       const hvl_float &gauss_base, const hvl_float &b1pb2, const hvl_float &reca3, const hvl_float &sqrtz)
{
	if (a3 == ctx->HVL_ZERO)
		return ToDouble(-ctx->HVL_HALF * a0 * Erf(sqrtz) * gauss_base);

	const hvl_float a3e = Exp(a3);
	const hvl_float nmr = gauss_base * a0 * (-reca3 + (a3e / (Pow(a3e - ctx->HVL_ONE, 2) * b1pb2)));

	return ToDouble(reca3 * nmr / b1pb2);
}

/* Oneshot HVL calculators with no precalculated parameters */
inline double HVL_t_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a0 = ToHVLFloat(a0);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);

	return HVL_t_internal(ctx, hvl_a0, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3));
}

inline double HVL_dx_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a0 = ToHVLFloat(a0);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);
	const hvl_float nmr = CalcNmr(ctx, hvl_a0, q, expSqrtz2);

	return HVL_dx_internal(ctx, hvl_a0, hvl_a2, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3), nmr, xmu);
}

inline double HVL_da0_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	(void)a0;

	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);

	return HVL_da0_internal(ctx, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3));
}

inline double HVL_da1_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a0 = ToHVLFloat(a0);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);
	const hvl_float nmr = CalcNmr(ctx, hvl_a0, q, expSqrtz2);

	return HVL_da1_internal(ctx, hvl_a0, hvl_a2, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3), nmr, xmu);
}

inline double HVL_da2_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a0 = ToHVLFloat(a0);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);
	const hvl_float nmr = CalcNmr(ctx, hvl_a0, q, expSqrtz2);

	return HVL_da2_internal(ctx, hvl_a0, hvl_a2, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3), nmr, xmu);
}

inline double HVL_da3_oneshot(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	const hvl_float hvl_x = ToHVLFloat(x);
	const hvl_float hvl_a0 = ToHVLFloat(a0);
	const hvl_float hvl_a1 = ToHVLFloat(a1);
	const hvl_float hvl_a2 = ToHVLFloat(a2);
	const hvl_float hvl_a3 = ToHVLFloat(a3);

	const hvl_float xmu = CalcXmu(hvl_x, hvl_a1);
	const hvl_float q = CalcQ(ctx, hvl_a2);
	const hvl_float sqrtz = CalcSqrtz(xmu, q);
	const hvl_float expSqrtz2 = CalcExpSqrtz2(sqrtz);
	const hvl_float b1pb2 = CalcB1pB2(ctx, hvl_a3, sqrtz);
	const hvl_float gauss_base = CalcGaussBase(ctx, q, expSqrtz2);

	return HVL_da3_internal(ctx, hvl_a0, hvl_a3, gauss_base, b1pb2, CalcReca3(ctx, hvl_a3), sqrtz);
}

typedef double (LIBHVL_DLLCALL *fpcalcfun)(const HVL_Context *, double, double, double, double, double); /*!< Pointer to function that does the calculation */
typedef struct {
	HVL_Pair *buffer;		/*!< Buffer where the results of the calculations will be stored */
	const HVL_Context *ctx;		/*!< Context of the current MPFR setup */
	double from;			/*!< Value of X to start counting from */
	double step;			/*!< Difference between two consecutive values of X */
	size_t iters;			/*!< Number of calculations to do */
	double a0;			/*!< HVL a0 parameter */
	double a1;			/*!< HVL a1 parameter */
	double a2;			/*!< HVL a2 parameter */
	double a3;			/*!< HVL a3 parameter */
	fpcalcfun calcfun;		/*!< Pointer to function that does the calculation */
} ThreadParams;

#ifdef LIBHVL_THREADING_WIN32
	#define HVL_THREAD HANDLE
#elif defined LIBHVL_THREADING_PTHREAD
	#define HVL_THREAD pthread_t
#endif

/** Multithreading helpers */

/*
 * @brief Allocates and initializes new HVLInternalValues structure
 *
 * @param[in] Number of x,y value pairs to store
 *
 * @return Pointer to struct HVLInternalValues
 */
static
HVL_Range * HVL_alloc_range(const size_t count)
{
	HVL_Range *pv;

	if (count == 0)
		return nullptr;

	try {
		pv = new HVL_Range;
	} catch (std::bad_alloc&) {
		return nullptr;
	}

	try {
		pv->p = new HVL_Pair[count];
	} catch (std::bad_alloc&) {
		delete pv;
		return nullptr;
	}
	pv->count = count;

	return pv;
}

static
size_t GetNumCPUs()
{
#ifdef LIBHVL_PLATFORM_WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#elif defined LIBHVL_PLATFORM_UNIX
	int cpus = get_nprocs();
	return cpus > 0 ? cpus : 1;
#endif //HVL_THREADING_
}

/*
 * @brief Returns the number of threads to use for calculation
 *
 * @return Number of threads to use for calculation
 */
static
size_t GetNumThreads(const size_t arraySize)
{
#ifdef LIBHVL_USE_CPP17
	const int cacheLineSize = std::hardware_destructive_interference_size;
#else
	const int cacheLineSize = 64; /* Assume 64-byte L1 cache line size */
#endif //LIBHVL_USE_CPP17
	const size_t numCPU = GetNumCPUs();
	const size_t maxChunks = arraySize * sizeof(double) / cacheLineSize;

	if (maxChunks < 1)
		return 1;
	if (maxChunks < numCPU)
		return maxChunks;
	return numCPU;
}

#ifdef LIBHVL_THREADING_PTHREAD
void ThreadCleanup(void *arg)
{
	(void)arg;
#ifdef LIBHVL_USE_MPFR
	mpfr_free_cache();
#endif // LIBHVL_USE_MPFR
}
#endif // LIBHVL_THREADING_PTHREAD

/*
 * @brief Worker thread function that calcualtes HVL function values for a given range
 *
 * @param arg Pointer to ThreadParams struct
 */
static
#ifdef LIBHVL_THREADING_WIN32
DWORD WINAPI WorkerFunc(void *arg)
#elif defined LIBHVL_THREADING_PTHREAD
void *WorkerFunc(void *arg)
#endif // LIBHVL_THREADING_
{
	const ThreadParams *tp = static_cast<const ThreadParams *>(arg);
	double x = tp->from;
	HVL_Pair *buf = tp->buffer;

	CheckPrecision(tp->ctx);
#ifdef LIBHVL_THREADING_PTHREAD
	pthread_cleanup_push(ThreadCleanup, nullptr);
#endif // LIBHVL_THREADING_PTHREAD

	for (size_t iter = 0; iter < tp->iters; iter++) {
		double r = tp->calcfun(tp->ctx, x, tp->a0, tp->a1, tp->a2, tp->a3);

		buf->x = x;
		buf->y = r;
		buf++;
		x += tp->step;
		#ifdef LIBHVL_THREADING_PTHREAD
		pthread_testcancel();
		#endif // LIBHVL_THREADING_PTHREAD
	}

#ifdef LIBHVL_USE_MPFR
	mpfr_free_cache();
#endif // LIBHVL_USE_MPFR

#ifdef LIBHVL_THREADING_PTHREAD
	pthread_cleanup_pop(0);
#endif // LIBHVL_THREADING_PTHREAD

	return 0;
}

/*
 * @brief Divides the HVL calculation between all available CPUs and waits for the calculation to finish
 *
 * @param[in,out] pv Double pointer to the HVLValues that are to be returned
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return True on success, false on failure
 */
static
HVL_RetCode LaunchWorkersAndWait(HVL_Range **pv, const HVL_Context *ctx, double from, double to, double step, double a0, double a1, double a2, double a3, fpcalcfun calcfun)
{
	HVL_THREAD *threads;
	ThreadParams *tps = nullptr;
	HVL_RetCode tRet;
	size_t idx;
	HVL_Pair *buffer;
#ifdef LIBHVL_THREADING_WIN32
	DWORD *threadIDs;
#elif defined LIBHVL_THREADING_PTHREAD
	int ret;
#else
	#error "No threading model has been specified"
#endif // LIBHVL_THREADING_
	const size_t count = (size_t)floor(((to - from) / step) + 0.5);
	if (count < 1)
		return HVL_E_INVALID_ARG;

	const size_t numThreads = GetNumThreads(count);
	const size_t itersPerThread = count / numThreads;

	try {
		threads = new HVL_THREAD[numThreads];
	} catch (std::bad_alloc&) {
		return HVL_E_NO_MEMORY;
	}
#ifdef LIBHVL_THREADING_WIN32
	try {
		threadIDs = new DWORD[numThreads];
	} catch (std::bad_alloc&) {
		delete[] threads;
		return HVL_E_NO_MEMORY;
	}
#endif // LIBHVL_THREADING_WIN32

	try {
		tps = new ThreadParams[numThreads];
	} catch (std::bad_alloc&) {
		tRet = HVL_E_NO_MEMORY;
		goto out;
	}

	*pv = HVL_alloc_range(count);
	if (*pv == nullptr) {
		tRet = HVL_E_NO_MEMORY;
		goto out;
	}

	buffer = (*pv)->p;

	for (idx = 0; idx < numThreads; idx++) {
		tps[idx].buffer = buffer + (itersPerThread * idx);
		tps[idx].ctx = ctx;
		tps[idx].from = from + (idx * step * itersPerThread);
		tps[idx].step = step;
		tps[idx].iters = itersPerThread;
		tps[idx].a0 = a0;
		tps[idx].a1 = a1;
		tps[idx].a2 = a2;
		tps[idx].a3 = a3;
		tps[idx].calcfun = calcfun;
	}
	/* Make an adjustment for the last thread */
	tps[numThreads - 1].iters = count - (itersPerThread * (numThreads - 1)); /* Make sure that we do not leave out anything due to rounding */

	MPFR_TRY_BEG
	/* Launch threads */
	for (idx = 0; idx < numThreads; idx++) {
	#ifdef LIBHVL_THREADING_WIN32
		threads[idx] = CreateThread(NULL, 0, WorkerFunc, &tps[idx], 0, &threadIDs[idx]);
		if (threads[idx] == 0)
			goto err_unwind;
	#elif defined LIBHVL_THREADING_PTHREAD
		ret = pthread_create(&threads[idx], NULL, WorkerFunc, &tps[idx]);
		if (ret)
			goto err_unwind;
	#endif // LIBHVL_THREADING_
	}

#ifdef LIBHVL_THREADING_WIN32
	{
		size_t thrIdx = 0;

		while (thrIdx < numThreads) {
			auto wRet = WaitForSingleObject(threads[thrIdx], 50);
			if (mpfr_assert_triggered)
				goto err_unwind;
			if (wRet == WAIT_OBJECT_0) {
				CloseHandle(threads[thrIdx]);
				thrIdx++;
			}
		}
	}
#elif defined LIBHVL_THREADING_PTHREAD
	#ifdef _WIN32
		for (size_t idx = 0; idx < numThreads; idx++)
			pthread_join(threads[idx], NULL);
	#else
	{

		struct timespec ts;
		ts.tv_sec = 0;
		ts.tv_nsec = 50000000;

		size_t thrIdx = 0;

		while (thrIdx < numThreads) {
			const auto jRet = pthread_timedjoin_np(threads[thrIdx], NULL, &ts);
			if (mpfr_assert_triggered)
				goto err_unwind;
			if (jRet == 0)
				thrIdx++;
		}
	}
	#endif // _WIN32
#endif // LIBHVL_THREADING_
	MPFR_TRY_END
	MPFR_CATCH_BEG
	goto err_unwind;
	MPFR_CATCH_END

	tRet = HVL_OK;
	goto out;

/* Something went wrong when we were starting the workers, terminate
   all workers that might have been started and exit */
err_unwind:
	if (mpfr_assert_triggered) {
		/* We are NOT reinstalling the original SIGABRT handler as some other threads
		 * that are still running may trigger the assert too and kill the whole program.
		 * In this case it should be the responsibility of the caller to reinstall
		 * the original SIGABRT handler should it need to.
		 * On the other hand a program that runs into this kind of MPFR exception
		 * should not try to do anything else but log what happened and terminate immediately.
		 * Safe recovery from this state is NOT POSSIBLE, do not even try!
		 */
		return HVL_E_MPFR_ASSERT;
	}

	tRet = HVL_E_INTERNAL;
	for (;;) {
	#ifdef LIBHVL_THREADING_WIN32
		TerminateThread(threads[idx - 1], WAIT_OBJECT_0);
		CloseHandle(threads[idx - 1]);
	#elif defined LIBHVL_THREADING_PTHREAD
		int _ret;
		void *retval;

		pthread_cancel(threads[idx - 1]);
		_ret = pthread_join(threads[idx - 1], &retval);
		if (_ret)
			abort();  /* Something has gone very wrong when killing the thread */
	#endif // LIBHVL_THREADING_
		if (idx == 1)
			break;
		else
			idx--;
	}
	HVL_free_range(*pv);
	*pv = nullptr;
out:
	delete[] tps;
	delete[] threads;
#ifdef LIBHVL_THREADING_WIN32
	delete[] threadIDs;
#endif // LIBHVL_THREADING_WIN32

	return tRet;
}

HVL_Context * MakeHVLContext(const int prec_bits)
{
	HVL_Context *ctx = nullptr;

	try {
		ctx = new HVL_Context(prec_bits);
	} catch (std::bad_alloc&) {
		ctx = nullptr;
	}

	return ctx;
}

HVL_Prepared * MakeHVLPrepared(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	HVL_Prepared *prep = nullptr;

	try {
		prep = new HVL_Prepared(ctx, x, a0, a1, a2, a3);
	} catch (std::bad_alloc&) {
		prep = nullptr;
	}

	return prep;
}

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/* Public library functions */

LIBHVL_DLLEXPORT void LIBHVL_DLLCALL
HVL_free_context(HVL_Context *ctx)
{
	delete ctx;
}

LIBHVL_DLLEXPORT void LIBHVL_DLLCALL
HVL_free_prepared(HVL_Prepared *prep)
{
	delete prep;
}

LIBHVL_DLLEXPORT void LIBHVL_DLLCALL
HVL_free_range(HVL_Range *r)
{
	if (r == nullptr)
		return;

	delete[] r->p;
	delete r;
}

LIBHVL_DLLEXPORT void LIBHVL_DLLCALL
HVL_free_threadlocal_cache()
{
#ifdef LIBHVL_USE_MPFR
	mpfr_free_cache();
#endif // LIBHVL_USE_MPFR
}

/*
 * @brief Creates HVL calculation context with the requested precision. Of an old context is passed as the input, it is free'd before a new one is set
 *
 * @param[in] ctx HVL calculation context
 * @param[in] prec Requested precision in significant digits
 *
 * @return HVL_OK on success, error code otherwise
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_make_context(HVL_Context **ctx, const int prec)
{
	#ifndef LIBHVL_USE_MPFR
		(void)prec;
	#endif // LIBHVL_USE_MPFR
	int prec_bits = 0;

	if (prec < 1)
		return HVL_E_INVALID_ARG;

	if (*ctx != nullptr)
		HVL_free_context(*ctx);

	#ifdef LIBHVL_USE_MPFR
		prec_bits = mpfr::digits2bits(prec);
		mpfr::mpreal::set_default_prec(prec_bits);
	#endif // LIBHVL_USE_MPFR

	*ctx = MakeHVLContext(prec_bits);
	if (*ctx == nullptr)
		return HVL_E_NO_MEMORY;

	return HVL_OK;
}

LIBHVL_DLLEXPORT enum HVL_RetCode LIBHVL_DLLCALL
HVL_prepare(HVL_Prepared **prep, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	*prep = MakeHVLPrepared(ctx, x, a0, a1, a2, a3);
	if (*prep == nullptr)
		return HVL_E_NO_MEMORY;

	HVL_Prepared *dprep = *prep;

	const_cast<hvl_float&>(dprep->xmu) = CalcXmu(dprep->x, dprep->a1);

	const hvl_float q = CalcQ(ctx, dprep->a2);
	const_cast<hvl_float&>(dprep->sqrtz) = CalcSqrtz(dprep->xmu, q);

	const hvl_float expSqrtz2 = CalcExpSqrtz2(dprep->sqrtz);
	const_cast<hvl_float&>(dprep->gauss_base) = CalcGaussBase(ctx, q, expSqrtz2);

	const_cast<hvl_float&>(dprep->b1pb2) = CalcB1pB2(ctx, dprep->a3, dprep->sqrtz);
	const_cast<hvl_float&>(dprep->nmr) = CalcNmr(ctx, dprep->a0, q, expSqrtz2);

	const_cast<hvl_float&>(dprep->reca3) = CalcReca3(ctx, dprep->a3);

	return HVL_OK;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_t_prepared(double *t, const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*t = HVL_t_internal(&prep->ctx, prep->a0, prep->a3, prep->gauss_base, prep->b1pb2, prep->reca3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_dx_prepared(double *dx, const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*dx = HVL_dx_internal(&prep->ctx, prep->a0, prep->a2, prep->a3,
				      prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da0_prepared(double *da0, const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da0 = HVL_da0_internal(&prep->ctx, prep->a3, prep->gauss_base, prep->b1pb2, prep->reca3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da1_prepared(double *da1,const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da1 = HVL_da1_internal(&prep->ctx,
				        prep->a0, prep->a2, prep->a3,
					prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da2_prepared(double *da2, const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da2 =  HVL_da2_internal(&prep->ctx,
					 prep->a0, prep->a2, prep->a3,
					 prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da3_prepared(double *da3, const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da3 =  HVL_da3_internal(&prep->ctx,
					 prep->a0, prep->a3,
					 prep->gauss_base, prep->b1pb2, prep->reca3, prep->sqrtz);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_t(double *t, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*t = HVL_t_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_dx(double *dx, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*dx = HVL_dx_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da0(double *da0, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da0 = HVL_da0_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da1(double *da1, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da1 = HVL_da1_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da2(double *da2, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da2 = HVL_da2_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}


LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da3(double *da3, const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	HVL_RetCode ret;
	MPFR_TRY_BEG
		*da3 = HVL_da3_oneshot(ctx, x, a0, a1, a2, a3);
		ret = HVL_OK;
	MPFR_TRY_END
	MPFR_CATCH_BEG
		ret = HVL_E_MPFR_ASSERT;
	MPFR_CATCH_END

	return ret;
}

/*
 * @brief Calculate values of HVL function for a given range of X
 *
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_t_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	    const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx, from, to, step, a0, a1, a2, a3, HVL_t_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

/*
* @brief Calculate values of HVL/dx function for a given range of X
*
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_dx_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	     const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx, from, to, step, a0, a1, a2, a3, HVL_dx_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

/*
 * @brief Calculate values of HVL/da0 function for a given range of X
 *
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da0_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx,  from, to, step, a0, a1, a2, a3, HVL_da0_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

/*
 * @brief Calculate values of HVL/da1 function for a given range of X
 *
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da1_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx, from, to, step, a0, a1, a2, a3, HVL_da1_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

/*
 * @brief Calculate values of HVL/da2 function for a given range of X
 *
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da2_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx, from, to, step, a0, a1, a2, a3, HVL_da2_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

/*
 * @brief Calculate values of HVL/da3 function for a given range of X
 *
 * @param[in,out] rng Range of computed values
 * @param[in] ctx Calculation context
 * @param[in] from Value of X to calculate from
 * @param[in] to Value of X to calculate to
 * @param[in] step Difference between two consecutive values of X
 * @param[in] a0 HVL a0 parameter
 * @param[in] a1 HVL a1 parameter
 * @param[in] a2 HVL a2 parameter
 * @param[in] a3 HVL a3 parameter
 *
 * @return Pointer to HVLInternalValues struct, NULL on failure
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_da3_range(HVL_Range **rng, const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	*rng = nullptr;
	if (to <= from)
		return HVL_E_INVALID_ARG;

	HVL_RetCode ret = LaunchWorkersAndWait(rng, ctx, from, to, step, a0, a1, a2, a3, HVL_da3_oneshot);
	if (ret != HVL_OK)
		HVL_free_range(*rng);
	return ret;
}

#ifdef __cplusplus
}
#endif // __cplusplus
