/*
 * conv.h
 *
 *  Created on: Oct 18, 2014
 *      Author: igkiou
 */

#ifndef CONV_H_
#define CONV_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "mkl_service.h"
#include "mkl_vsl.h"

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
