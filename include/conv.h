/*
 * conv.h
 *
 *  Created on: Oct 18, 2014
 *      Author: igkiou
 */

#ifndef CONV_H_
#define CONV_H_

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

void dconv1();

void dconv2();

void dconvn();

//void seed(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info);

void drunif(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info);

void drnorm(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info);

void rngDestroy(RngEngine* rng, RngErrorType* info);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* CONV_H_ */
