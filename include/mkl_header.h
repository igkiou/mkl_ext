/*
 * mkl_ext//mkl_ext/include/use_mkl.h/use_mkl.h
 *
 *  Created on: Jun 6, 2013
 *      Author: igkiou
 */

#ifndef MKL_HEADER_H_
#define MKL_HEADER_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <mkl.h>
#include <omp.h>

typedef MKL_INT BlasInt;
#define BLASFUNC(NAME) NAME

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MKL_HEADER_H_ */
