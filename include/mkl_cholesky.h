/*
 * cholesky.h
 *
 *  Created on: Mar 15, 2012
 *      Author: igkiou
 */

#ifndef __CHOLESKY_H__
#define __CHOLESKY_H__

#include <mkl.h>
#include <omp.h>

void dchud(char *uplo, MKL_INT *pN, double *a, double *work, MKL_INT *info);

void dchdd(char *uplo, MKL_INT *pN, double *a, double *work, MKL_INT *info);

void dchr(char *uplo, MKL_INT *N, double *alpha, double *x, MKL_INT *incx, \
		double *a, MKL_INT *lda, double *work, MKL_INT *lwork, MKL_INT *info);

void dchmv(char *uplo, MKL_INT *N, double *a, MKL_INT *lda, double *x, \
			MKL_INT *incx);

void dchrk(char *uplo, char *trans, MKL_INT *N, MKL_INT *K, double *alpha, \
			double *a, MKL_INT *lda, double *c, MKL_INT *ldc, double *work, \
			MKL_INT *lwork, MKL_INT *info);


void dchmm(char *uplo, char *side, MKL_INT *M, MKL_INT *N, double *alpha, \
			double *a, MKL_INT *lda, double *b, MKL_INT *ldb);

//void CHOA();
//void CHOD();
//void CHOAK();
//void CHODK();

#endif /* __CHOLESKY_H__ */
