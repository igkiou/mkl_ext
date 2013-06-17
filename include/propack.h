/*
 * mkl_ext//mkl_ext/include/propack.h/propack.h
 *
 *  Created on: Jun 7, 2013
 *      Author: igkiou
 */

#ifndef PROPACK_H_
#define PROPACK_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include "blas_header.h"

/*
 * TODO: Maybe add lda and lwork arguments for uniformity with BLAS calls? In
 * that case, also add argument checks and errors in info.
 */

void dlansvd(char* jobu, char* jobvt, BlasInt* m, BlasInt* n, double* a,
			double* u, double* vt, BlasInt* k, BlasInt* job);

void dlanbpro()

//void dchdk(double* a, BlasInt* n, BlasInt* k, BlasInt* m, double* work);
//void dcha(double* a, BlasInt* n, BlasInt* k, double* work);
//void dchak(double* a, BlasInt* n, BlasInt* k, BlasInt* m, double* work);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* PROPACK_H_ */
