#include "hvl.hpp"

#if ECHMET_MATH_HVL_DLLBUILD == 0

// -> EOF

#else

#include <limits>

#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	#include <algorithm>
#endif
#include "mpreal.h"

//-----------------------------------------------------------------------------
#if ECHMET_MATH_HVL_DLLBUILD == -1
	#define DLL_EXPORT
#else
#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	#define DLL_EXPORT __declspec(dllexport)
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	#define DLL_EXPORT __attribute__ ((visibility("default")))
#endif // HVL_PLATFORM_

#endif
	
//-----------------------------------------------------------------------------
#define ECHMET_MATH_HVL_TODOUBLE(VAL) (VAL).toDouble()
#define Infinity mpfr::const_infinity()
typedef mpfr::mpreal float_type;

struct HVLContext {
	float_type HVL_SQRT2;
	float_type HVL_PI;
	float_type HVL_SQRTPI;
	float_type HVL_05;
	float_type HVL_00;
	float_type HVL_1;
	float_type HVL_2;
	float_type HVL_SQRT2REC;
};

//=============================================================================
// CODE

//---------------------- Internal calculation functions -----------------------

inline float_type Ln(float_type x)
{
	return mpfr::log(x);
}

//-----------------------------------------------------------------------------
inline float_type Exp(float_type x)
{
	static const float_type maxx = Ln( std::numeric_limits<float_type>::max() );

	return x > maxx ? Infinity : mpfr::exp(x);
}

//-----------------------------------------------------------------------------
inline float_type Erf(float_type x)
{
	return mpfr::erf(x);
}

//-----------------------------------------------------------------------------
inline float_type Pow(float_type x, float_type p)
{
	return mpfr::pow(x, p);
}

//-----------------------------------------------------------------------------
inline float_type Sqrt(float_type x)
{
	return mpfr::sqrt(x);
}

//-----------------------------------------------------------------------------
// DLL Code

	namespace echmet {
		namespace math {
#if ECHMET_MATH_HVL_DLLBUILD == -1
	/*
	 * When building as a module nest everything into
	 * echmet::math::impl namespace
	*/
		/*
		 * Include internal types when the code is
		 * not built as a library
		 */
		#include "hvl_types.hpp"
			namespace impl {
#elif ECHMET_MATH_HVL_DLLBUILD == 1
	/* Make sure that we export all functions with C linkage */
	#ifdef __cplusplus
	extern "C" {
	#endif //__cplusplus
#endif //ECHMET_MATH_HVL_DLLBUILD

/*
 * @brief Calculate value of HVL function
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVL(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{
	float_type result;
	float_type gauss;

	const float_type u      = a1;
	const float_type s      = a2;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	gauss = a0 * q / ctx->HVL_SQRTPI * Exp(-sqrtz*sqrtz);

	if (a3 == 0) return ECHMET_MATH_HVL_TODOUBLE(gauss);

	// ! Conditional return above

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	result = ctx->HVL_1 / a3 * gauss / (b1 + b2);

	return ECHMET_MATH_HVL_TODOUBLE(result);

}

//-----------------------------------------------------------------------------

/*
 * @brief Calculate value of HVL/dx
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL/dx function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVLdx(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{

	float_type result;
	float_type gauss_der;

	const float_type u      = a1;
	const float_type s      = a2;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	gauss_der = -a0 * (x - u) * q / ( s * s ) / ctx->HVL_SQRTPI * Exp(-sqrtz*sqrtz); // Gauss derivative

	if (a3 == 0) return ECHMET_MATH_HVL_TODOUBLE(gauss_der);

	// ! Conditional return above

	const float_type nmr    = -a0 * q * q * Exp(-sqrtz*sqrtz) * Exp(-sqrtz*sqrtz) / ctx->HVL_PI;

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	result = ( ctx->HVL_1 / a3 ) * ( ( gauss_der / (b1 + b2) ) + ( nmr / ( (b1 + b2) * (b1 + b2) ) ) );

	return ECHMET_MATH_HVL_TODOUBLE(result);

}

//-----------------------------------------------------------------------------
/*
 * @brief Calculate value of HVL/da0
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL/da0 function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVLda0(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{

	float_type result;
	float_type gauss_der;

	const float_type u      = a1;
	const float_type s      = a2;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	gauss_der = q * Exp(-sqrtz*sqrtz) / ctx->HVL_SQRTPI; // Gauss derivative

	if (a3 == 0) return ECHMET_MATH_HVL_TODOUBLE(gauss_der);

	// ! Conditional return above

	//const float_type nmr    = -a0 * q / s / HVL_PI * Exp(-sqrtz*sqrtz) * Exp(-sqrtz*sqrtz);

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	result = ( ctx->HVL_1 / a3 ) * ( gauss_der / (b1 + b2) );

	return ECHMET_MATH_HVL_TODOUBLE(result);

}

//-----------------------------------------------------------------------------
/*
 * @brief Calculate value of HVL/da1
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL/da1 function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVLda1(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{

	float_type result;
	float_type gauss_der;

	const float_type u      = a1;
	const float_type s      = a2;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	gauss_der = a0 * (x - u) * q * Exp(-sqrtz*sqrtz) / ctx->HVL_SQRTPI /s /s; // Gauss derivative

	if (a3 == 0) return ECHMET_MATH_HVL_TODOUBLE(gauss_der);

	// ! Conditional return above

	const float_type nmr    = a0 * q * q * Exp(-sqrtz*sqrtz) * Exp(-sqrtz*sqrtz) / ctx->HVL_PI;

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	result = ( ctx->HVL_1 / a3 ) * ( ( gauss_der / (b1 + b2) ) + ( nmr / ( (b1 + b2) * (b1 + b2) ) ) );

	return ECHMET_MATH_HVL_TODOUBLE(result);

}

//-----------------------------------------------------------------------------
/*
 * @brief Calculate value of HVL/da2
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL/da2 function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVLda2(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{

	float_type result;
	float_type gauss_der;

	const float_type u      = a1;
	const float_type s      = a2;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	gauss_der = a0 * ( (x - u) * (x - u) * q * Exp(-sqrtz*sqrtz) / ( s * s * s ) / ctx->HVL_SQRTPI - ( q * Exp(-sqrtz*sqrtz) / s / ctx->HVL_SQRTPI ) );
	// Gauss derivative

	if (a3 == 0) return ECHMET_MATH_HVL_TODOUBLE(gauss_der);

	// ! Conditional return above

	const float_type nmr    = a0 * q * q * Exp(-sqrtz*sqrtz) * Exp(-sqrtz*sqrtz) * (x - u) / ctx->HVL_PI /s;

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	result = ( ctx->HVL_1 / a3 ) * ( ( gauss_der / (b1 + b2) ) + ( nmr / ( (b1 + b2) * (b1 + b2) ) ) );

	return ECHMET_MATH_HVL_TODOUBLE(result);
}

//-----------------------------------------------------------------------------
/*
 * @brief Calculate value of HVL/da3
 *
 * @param[in] ctx Calculation context
 * @param[in] x Independent variable
 * @param[in] a0 Parameter a0
 * @param[in] a1 Parameter a1
 * @param[in] a2 Parameter a2
 * @param[in] a3 Parameter a3
 *
 * @return Value of HVL/da3 function
 */
DLL_EXPORT double ECHMET_MATH_HVL_DLLCALL
HVLda3(const HVLContext *ctx, double x, double a0, double a1, double a2, double a3)
{

	float_type result;

	const float_type u      = a1;
	const float_type s      = a2;
	const float_type a      = a3;

	const float_type q      = ctx->HVL_SQRT2REC / s;
	const float_type sqrtz  = (x - u) * q;

	if (a3 == 0)  {

		result = -ctx->HVL_05 * a0 * Erf(sqrtz) * q * Exp(-sqrtz*sqrtz) / ctx->HVL_SQRTPI;
		return ECHMET_MATH_HVL_TODOUBLE(result);

	}

	const float_type b1     = ctx->HVL_1 / ( Exp(a3) - ctx->HVL_1);
	const float_type b2     = ctx->HVL_05 * ( ctx->HVL_1 + Erf(sqrtz) );

	const float_type nmr    = a0 * q * Exp(-sqrtz*sqrtz) * ( -1 / a + ( Exp(a) / ( (Exp(a) - 1) * (Exp(a) - 1) ) / ( b1 + b2 ) ) ) / ctx->HVL_SQRTPI;

	result = ( ctx->HVL_1 / a ) * nmr / ( b1 + b2 );

	return ECHMET_MATH_HVL_TODOUBLE(result);

}

typedef double (ECHMET_MATH_HVL_DLLCALL *fpcalcfun)(const HVLContext *, double, double, double, double, double); /*!< Pointer to function that does the calculation */
typedef struct {
	HVLInternalPair *buffer;	/*!< Buffer where the results of the calculations will be stored */
	const HVLContext *ctx; /*!< Context of the current MPFR setup */
	double from;	/*!< Value of X to start counting from */
	double step;	/*!< Difference between two consecutive values of X */
	size_t iters;	/*!< Number of calculations to do */
	double a0;		/*!< HVL a0 parameter */
	double a1;		/*!< HVL a1 parameter */
	double a2;		/*!< HVL a2 parameter */
	double a3;		/*!< HVL a3 parameter */
	fpcalcfun calcfun; /*!< Pointer to function that does the calculation */
} ThreadParams;

#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	#define HVL_THREAD HANDLE
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	#define HVL_THREAD pthread_t
#endif

/*
 * @brief Returns the number of threads to use for calculation
 *
 * @return Number of threads to use for calculation
 */
size_t GetNumThreads()
{
#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	int cpus = get_nprocs();
	return cpus > 0 ? cpus : 1;
#endif //HVL_PLATFORM_
}

/*
 * @brief Allocates and initializes new HVLInternalValues structure
 *
 * @param[in] Number of x,y value pairs to store
 *
 * @return Pointer to struct HVLInternalValues
 */
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLAlloc(const size_t count)
{
	HVLInternalValues *pv;

	if (count == 0)
		return NULL;

	pv = (HVLInternalValues *)malloc(sizeof(HVLInternalValues));
	if (pv == NULL)
		return NULL;

	pv->p = (HVLInternalPair *)malloc(sizeof(HVLInternalPair) * count);
	if (pv->p == NULL) {
		free(pv);
		return NULL;
	}
	pv->count = count;

	return pv;
}

/*
 * @brief Worker thread function that calcualtes HVL function values for a given range
 *
 * @param arg Pointer to ThreadParams struct
 */
#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
DWORD WINAPI WorkerFunc(void *arg)
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
void *WorkerFunc(void *arg)
#endif
{
	ThreadParams *tp = (ThreadParams *)arg;
	double x = tp->from;
	HVLInternalPair *buf = tp->buffer;
#ifdef ECHMET_MATH_HVL_PLATFORM_UNIX
	int old;
	pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &old);
#endif

	for (size_t iter = 0; iter < tp->iters; iter++) {
		double r = tp->calcfun(tp->ctx, x, tp->a0, tp->a1, tp->a2, tp->a3);

		buf->x = x;
		buf->y = r;
		buf++;
		x += tp->step;
	}

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
bool ECHMET_MATH_HVL_DLLCALL LaunchWorkersAndWait(HVLInternalValues **pv, const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3, fpcalcfun calcfun)
{
	HVL_THREAD *threads;
	ThreadParams *tps;
	bool bret;
	size_t idx;
#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	DWORD *threadIDs;
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	int ret;
#endif
	size_t numCPU = GetNumThreads();
	size_t count = (size_t)floor(((to - from) / step) + 0.5);
	size_t itersPerThread = count / numCPU;

	threads = (HVL_THREAD *)malloc(sizeof(HVL_THREAD) * numCPU);
	if (threads == NULL)
		return false;
#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	threadIDs = (DWORD *)malloc(sizeof(DWORD) * numCPU);
	if (threadIDs == NULL) {
		free(threads);
		return false;
	}
#endif

	tps = (ThreadParams *)malloc(sizeof(ThreadParams) * numCPU);
	if (tps == NULL) {
		bret = false;
		goto out;
	}

	*pv = HVLAlloc(count);
	if (*pv == NULL) {
		bret = false;
		goto out;
	}

	for (idx = 0; idx < numCPU-1; idx++) {
		tps[idx].buffer = &((*pv)->p[(itersPerThread * idx)]);
		tps[idx].ctx = ctx;
		tps[idx].from = from + (idx * step * itersPerThread);
		tps[idx].step = step;
		tps[idx].iters = itersPerThread;
		tps[idx].a0 = a0;
		tps[idx].a1 = a1;
		tps[idx].a2 = a2;
		tps[idx].a3 = a3;
		tps[idx].calcfun = calcfun;

	#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
		threads[idx] = CreateThread(NULL, 0, WorkerFunc, &tps[idx], 0, &threadIDs[idx]);
		if (threads[idx] == 0)
			goto err_unwind;
	#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
		ret = pthread_create(&threads[idx], NULL, WorkerFunc, &tps[idx]);
		if (ret)
			goto err_unwind;
	#endif
	}
	/* Last thread */
	tps[idx].buffer = &((*pv)->p[idx * itersPerThread]);
	tps[idx].ctx = ctx;
	tps[idx].from = from + (idx * step * itersPerThread);
	tps[idx].step = step;
	tps[idx].iters = count - (itersPerThread * idx); /* Make sure that we do not leave out anything due to rounding */
	tps[idx].a0 = a0;
	tps[idx].a1 = a1;
	tps[idx].a2 = a2;
	tps[idx].a3 = a3;
	tps[idx].calcfun = calcfun;

#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	threads[idx] = CreateThread(NULL, 0, WorkerFunc, &tps[idx], 0, &threadIDs[idx]);
	if (threads[idx] == 0)
		goto err_unwind;
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	ret = pthread_create(&threads[idx], NULL, WorkerFunc, &tps[idx]);
	if (ret)
		goto err_unwind;
#endif

#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
	WaitForMultipleObjects(numCPU, threads, TRUE, INFINITE);
	for (size_t idx = 0; idx < numCPU; idx++)
		CloseHandle(threads[idx]);
#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
	for (size_t idx = 0; idx < numCPU; idx++)
		pthread_join(threads[idx], NULL);
#endif
	bret = true;
	goto out;

/* Something went wrong when we were starting the workers, terminate
   all workers that might have been started and exit */
err_unwind:
	bret = false;
	for (;;) {
	#ifdef ECHMET_MATH_HVL_PLATFORM_WIN
		TerminateThread(threads[idx], WAIT_OBJECT_0);
		CloseHandle(threads[idx]);
	#elif defined ECHMET_MATH_HVL_PLATFORM_UNIX
		int _ret;
		void *retval;

		pthread_cancel(threads[idx]);
		_ret = pthread_join(threads[idx], &retval);
		if (_ret || (retval != PTHREAD_CANCELED))
			abort();  /* Something has gone very wrong when killing the thread */
	#endif
		if (idx == 0)
			break;
		else
			idx--;
	}
	free(*pv);
	*pv = NULL;
out:
	free(tps);
	free(threads);
#if ECHMET_MATH_HVL_PLATFORM_WIN
	free(threadIDs);
#endif

	return bret;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVL_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{
	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx, from, to, step, a0, a1, a2, a3, HVL) == false)
		return NULL;
	return values;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLdx_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx, from, to, step, a0, a1, a2, a3, HVLdx) == false)
		return NULL;
	return values;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLda0_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx,  from, to, step, a0, a1, a2, a3, HVLda0) == false)
		return NULL;
	return values;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLda1_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{
	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx, from, to, step, a0, a1, a2, a3, HVLda1) == false)
		return NULL;
	return values;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLda2_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{
	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx, from, to, step, a0, a1, a2, a3, HVLda2) == false)
		return NULL;
	return values;
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
DLL_EXPORT HVLInternalValues * ECHMET_MATH_HVL_DLLCALL
HVLda3_range(const HVLContext *ctx, double from, double to, double step, double a0, double a1, double a2, double a3)
{
	HVLInternalValues *values = NULL;
	if (to <= from)
		return NULL;

	if (LaunchWorkersAndWait(&values, ctx, from, to, step, a0, a1, a2, a3, HVLda3) == false)
		return NULL;
	return values;
}

/*
 * @brief Correctly frees HVLInternalValues struct
 *
 * @param[in] ptr Pointer to the HVLInternalValues struct
 */
DLL_EXPORT void ECHMET_MATH_HVL_DLLCALL
HVLFree(HVLInternalValues *ptr)
{
	free(ptr->p);
	free(ptr);
}

/*
 * @brief Sets MPFR library precision and creates calculation context if necessary
 *
 * @param[in] ctx Calculation context
 * @param[in] prec Precision in significant digits
 *
 * @return True on success, false if the calculation context could not have
 *	   been allocated
 */
DLL_EXPORT bool ECHMET_MATH_HVL_DLLCALL
HVLSetPrec(HVLContext **ctx, const int prec)
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(prec));

	if (*ctx == NULL) {
		try {
			*ctx = new HVLContext();
		} catch (std::bad_alloc& ) {
			return false;
		}
	}

	(*ctx)->HVL_SQRT2 = mpfr::sqrt(mpfr::mpreal("2.0"));
	(*ctx)->HVL_PI = mpfr::const_pi();
	(*ctx)->HVL_SQRTPI = mpfr::sqrt((*ctx)->HVL_PI);
	(*ctx)->HVL_05 = mpfr::mpreal("0.5");
	(*ctx)->HVL_00 = mpfr::mpreal("0.0");
	(*ctx)->HVL_1 = mpfr::mpreal("1.0");
	(*ctx)->HVL_2 = mpfr::mpreal(2);
	(*ctx)->HVL_SQRT2REC = (*ctx)->HVL_1 / (*ctx)->HVL_SQRT2;

	return true;
}

/*
 * @brief Frees calculation context
 *
 * @param[in] ctx Calculation context to be free'd
 */
DLL_EXPORT void ECHMET_MATH_HVL_DLLCALL
HVLFreeContext(HVLContext *ctx)
{
	if (ctx == NULL)
		return;
	delete ctx;
}

//-----------------------------------------------------------------------------

#if ECHMET_MATH_HVL_DLLBUILD == -1

			} //impl

#elif ECHMET_MATH_HVL_DLLBUILD == 1
	#ifdef __cplusplus
	}
	#endif //_cplusplus
#endif // ECHMET_MATH_HVL_DLLBUILD
		} //math
	} // echmet

#endif // ECHMET_MATH_HVL_DLLBUILD == 0
