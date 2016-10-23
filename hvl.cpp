#include "hvl.h"
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
#ifdef LIBHVL_THREADING_WIN32
DWORD WINAPI WorkerFunc(void *arg)
#elif defined LIBHVL_THREADING_PTHREAD
void *WorkerFunc(void *arg)
#endif // LIBHVL_THREADING_
{
	ThreadParams *tp = (ThreadParams *)arg;
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
bool LaunchWorkersAndWait(HVL_Range **pv, const HVL_Context *ctx, double from, double to, double step, double a0, double a1, double a2, double a3, fpcalcfun calcfun)
{
	HVL_THREAD *threads;
	ThreadParams *tps = nullptr;
	bool bret;
	size_t idx;
	HVL_Pair *buffer;
#ifdef LIBHVL_THREADING_WIN32
	DWORD *threadIDs;
#elif defined LIBHVL_THREADING_PTHREAD
	int ret;
#else
	#error "No threading model has been specified"
#endif // LIBHVL_THREADING_
	size_t count = (size_t)floor(((to - from) / step) + 0.5);
	size_t numThreads = GetNumThreads(count);
	size_t itersPerThread = count / numThreads;

	try {
		threads = new HVL_THREAD[numThreads];
	} catch (std::bad_alloc&) {
		return false;
	}
#ifdef LIBHVL_THREADING_WIN32
	try {
		threadIDs = new DWORD[numThreads];
	} catch (std::bad_alloc&) {
		delete[] threads;
		return false;
	}
#endif // LIBHVL_THREADING_WIN32

	try {
		tps = new ThreadParams[numThreads];
	} catch (std::bad_alloc&) {
		bret = false;
		goto out;
	}

	*pv = HVL_alloc_range(count);
	if (*pv == nullptr) {
		bret = false;
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
	WaitForMultipleObjects(numThreads, threads, TRUE, INFINITE);
	for (idx = 0; idx < numThreads; idx++)
		CloseHandle(threads[idx]);
#elif defined LIBHVL_THREADING_PTHREAD
	for (size_t idx = 0; idx < numThreads; idx++)
		pthread_join(threads[idx], NULL);
#endif // LIBHVL_THREADING_
	bret = true;
	goto out;

/* Something went wrong when we were starting the workers, terminate
   all workers that might have been started and exit */
err_unwind:
	bret = false;
	for (;;) {
	#ifdef LIBHVL_THREADING_WIN32
		TerminateThread(threads[idx], WAIT_OBJECT_0);
		CloseHandle(threads[idx]);
	#elif defined LIBHVL_THREADING_PTHREAD
		int _ret;
		void *retval;

		pthread_cancel(threads[idx]);
		_ret = pthread_join(threads[idx], &retval);
		if (_ret)
			abort();  /* Something has gone very wrong when killing the thread */
	#endif // LIBHVL_THREADING_
		if (idx == 0)
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

	return bret;
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

/*
 * @brief Creates HVL calculation context with the requested precision. Of an old context is passed as the input, it is free'd before a new one is set
 *
 * @param[in] ctx HVL calculation context
 * @param[in] prec Requested precision in significant digits
 *
 * @return HVL_OK on success, error code otherwise
 */
LIBHVL_DLLEXPORT HVL_RetCode LIBHVL_DLLCALL
HVL_make_context (HVL_Context **ctx, const int prec)
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

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_t_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_t_internal(&prep->ctx, prep->a0, prep->a3, prep->gauss_base, prep->b1pb2, prep->reca3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_dx_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_dx_internal(&prep->ctx, prep->a0, prep->a2, prep->a3,
			       prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da0_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_da0_internal(&prep->ctx, prep->a3, prep->gauss_base, prep->b1pb2, prep->reca3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da1_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_da1_internal(&prep->ctx,
			        prep->a0, prep->a2, prep->a3,
				prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da2_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_da2_internal(&prep->ctx,
				prep->a0, prep->a2, prep->a3,
				prep->gauss_base, prep->b1pb2, prep->reca3, prep->nmr, prep->xmu);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da3_prepared(const HVL_Prepared *prep)
{
	CheckPrecision(&prep->ctx);
	return HVL_da3_internal(&prep->ctx,
				prep->a0, prep->a3,
				prep->gauss_base, prep->b1pb2, prep->reca3, prep->sqrtz);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_t(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	return HVL_t_oneshot(ctx, x, a0, a1, a2, a3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_dx(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	return HVL_dx_oneshot(ctx, x, a0, a1, a2, a3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da0(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	return HVL_da0_oneshot(ctx, x, a0, a1, a2, a3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da1(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	return HVL_da1_oneshot(ctx, x, a0, a1, a2, a3);
}

LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da2(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	CheckPrecision(ctx);

	return HVL_da2_oneshot(ctx, x, a0, a1, a2, a3);
}


LIBHVL_DLLEXPORT double LIBHVL_DLLCALL
HVL_da3(const HVL_Context *ctx, const double x, const double a0, const double a1, const double a2, const double a3)
{
	return HVL_da3_oneshot(ctx, x, a0, a1, a2, a3);
}

/*
 * @brief Calculate values of HVL function for a given range of X
 *
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_t_range(const HVL_Context *ctx, const double from, const double to, const double step,
	    const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx, from, to, step, a0, a1, a2, a3, HVL_t_oneshot) == false)
		return nullptr;
	return r;
}

/*
* @brief Calculate values of HVL/dx function for a given range of X
*
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_dx_range(const HVL_Context *ctx, const double from, const double to, const double step,
	     const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx, from, to, step, a0, a1, a2, a3, HVL_dx_oneshot) == false)
		return nullptr;
	return r;
}

/*
 * @brief Calculate values of HVL/da0 function for a given range of X
 *
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_da0_range(const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx,  from, to, step, a0, a1, a2, a3, HVL_da0_oneshot) == false)
		return nullptr;
	return r;
}

/*
 * @brief Calculate values of HVL/da1 function for a given range of X
 *
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_da1_range(const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx, from, to, step, a0, a1, a2, a3, HVL_da1_oneshot) == false)
		return nullptr;
	return r;
}

/*
 * @brief Calculate values of HVL/da2 function for a given range of X
 *
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_da2_range(const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx, from, to, step, a0, a1, a2, a3, HVL_da2_oneshot) == false)
		return nullptr;
	return r;
}

/*
 * @brief Calculate values of HVL/da3 function for a given range of X
 *
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
LIBHVL_DLLEXPORT HVL_Range * LIBHVL_DLLCALL
HVL_da3_range(const HVL_Context *ctx, const double from, const double to, const double step,
	      const double a0, const double a1, const double a2, const double a3)
{
	HVL_Range *r = nullptr;
	if (to <= from)
		return nullptr;

	if (LaunchWorkersAndWait(&r, ctx, from, to, step, a0, a1, a2, a3, HVL_da3_oneshot) == false)
		return nullptr;
	return r;
}

#ifdef __cplusplus
}
#endif // __cplusplus
