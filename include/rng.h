/*
 * mkl_ext//mkl_ext/include/rng_header.h/rng_header.h
 *
 *  Created on: Jun 17, 2013
 *      Author: igkiou
 */

#ifndef RNG_H_
#define RNG_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if !defined(USE_RNG_VSL) && !defined(USE_RNG_SFMT) && !defined(USE_RNG_DSFMT)
#define USE_RNG_VSL
#endif
#if !defined(USE_RNG_MARSAGLIA) && !defined(USE_RNG_BOX_MULLER)
#define USE_RNG_BOX_MULLER
#endif

#if defined(USE_RNG_SFMT)
#define SFMT_MEXP 19937
#include "SFMT.h"
typedef uint32_t RngSeedType;
typedef int RngErrorType;
typedef struct RngEngine {
	sfmt_t m_sfmt;
} RngEngine;
#elif defined(USE_RNG_DSFMT)
#define DSFMT_MEXP 19937
#include "dSFMT.h"
typedef uint32_t RngSeedType;
typedef int RngErrorType;
typedef struct RngEngine {
	dsfmt_t m_dsfmt;
} RngEngine;
#elif defined(USE_RNG_VSL)
#include "mkl_service.h"
#include "mkl_vsl.h"
const MKL_INT kVSLBRNGMethod = VSL_BRNG_SFMT19937;
typedef MKL_INT RngSeedType;
typedef int RngErrorType;
typedef struct RngEngine {
	VSLStreamStatePtr m_stream;
} RngEngine;
#endif

#include "blas_header.h"

void rngInit(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info);

//void seed(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info);

void drunif(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info);

void drnorm(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info);

void rngDestroy(RngEngine* rng, RngErrorType* info);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* RNG_H_ */
