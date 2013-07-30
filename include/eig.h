/*
 * mkl_ext//mkl_ext/include/mkl_chol.h/eig.h
 *
 *  Created on: July 23, 2013
 *      Author: igkiou
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include "blas_header.h"

/*
 * TODO: Maybe add lda and lwork arguments for uniformity with BLAS calls? In
 * that case, also add argument checks and errors in info.
 */

void dchud_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info);

void dchud(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info);

void dchdd_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info);

void dchdd(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info);

void dchr_sub(char* uplo, BlasInt* n, double* alpha, double* a, double* work,
			BlasInt* info);

void dchr(char* uplo, BlasInt* n, double* alpha, double* x, BlasInt* incx,
		double* a, double* work, BlasInt* info);

void dchmv(char* uplo, BlasInt* n, double* a, double* x, BlasInt* incx);

void dchrk(char* uplo, char* trans, BlasInt* n, BlasInt* k, double* alpha,
		double* a, double* c, double* work, BlasInt* info);

void dchmm(char* side, char* uplo, BlasInt* m, BlasInt* n, double* alpha,
		double* a, double* b);

void dchex(double* a, BlasInt* n, BlasInt* k, BlasInt* l, double* work,
		BlasInt* job);

void dchd(double* a, BlasInt* n, BlasInt* k, double* work);

//void dchdk(double* a, BlasInt* n, BlasInt* k, BlasInt* m, double* work);
//void dcha(double* a, BlasInt* n, BlasInt* k, double* work);
//void dchak(double* a, BlasInt* n, BlasInt* k, BlasInt* m, double* work);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* CHOLESKY_H_*/
