/*
 * mutils.h
 *
 *  Created on: Mar 21, 2012
 *      Author: igkiou
 */

#ifndef UTILS_H_
#define UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "blas_header.h"

void ddimm(char* side, char* trans, BlasInt* m, BlasInt* n, double* alpha,
		double* a, double* b, double* beta, double* c);

void ddiag(char* trans, BlasInt* m, BlasInt* n, double* a, double* b);

double dtrac(BlasInt* m, BlasInt* n, double* a);

double dvsum(BlasInt* n, double* x, BlasInt* incx);

double dvprd(BlasInt* n, double* x, BlasInt* incx);

void dgeva(char* trans, BlasInt* n, BlasInt* m, double* alpha, double* a,
		double* b);

void dgevm(char* trans, BlasInt* m, BlasInt* n, double* alpha, double* a,
		double* b);

void dgesu(char* trans, BlasInt* n, BlasInt* m, double* alpha, double* a,
		double* beta, double* b);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* UTILS_H_ */
