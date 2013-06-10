/*
 * mkl_ext//mkl_ext/include/mkl_chol.h/mkl_chol.h
 *
 *  Created on: May 17, 2013
 *      Author: igkiou
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include "blas_header.h"

void dchud_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info);

void dchud(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info);

void dchdd_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info);

void dchdd(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info);

//void dchr(char* uplo, BlasInt* n, double* alpha, double* x, BlasInt* incx,
//		double* a, BlasInt* lda, double* work, BlasInt* lwork, BlasInt* info);
//
//void dchmv(char* uplo, BlasInt* n, double* a, BlasInt* lda, double* x,
//		BlasInt* incx);
//
//void dchrk(char* uplo, char* trans, BlasInt* n, BlasInt* K, double* alpha,
//		double* a, BlasInt* lda, double* c, BlasInt* ldc, double* work,
//		BlasInt* lwork, BlasInt* info);
//
//void dchmm(char* uplo, char* side, BlasInt* m, BlasInt* n, double* alpha,
//		double* a, BlasInt* lda, double* b, BlasInt* ldb);

//void CHOA();
//void CHOD();
//void CHOAK();
//void CHODK();

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* CHOLESKY_H_*/
