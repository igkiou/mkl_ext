/*
 * mkl_ext//mkl_ext/include/blas_header.h/blas_header.h
 *
 *  Created on: Jun 7, 2013
 *      Author: igkiou
 */

#ifndef BLAS_HEADER_H_
#define BLAS_HEADER_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef USE_BLAS_MKL
#include <mkl.h>
#include <omp.h>

typedef MKL_INT BlasInt;
#define BLASFUNC(NAME) NAME
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* BLAS_HEADER_H_ */
