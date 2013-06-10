/*
 * utils.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: igkiou
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

inline void datacpy(double* target, double* source, BlasInt size) {
	BlasInt incx = 1;
	dcopy(&size, source, &incx, target, &incx);
}

void ddimm(char* side, char* trans, BlasInt* m, BlasInt* n, double* alpha,
		double* a, double* b, double* beta, double* c) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt incx;
	BlasInt iter;
	double alphaval;
	if (*beta == 0) {
		BlasInt SCALN;
		datacpy(c, b, M * N);
		if (*side == 'L') {
			SCALN = N;
			incx = M;
			if (*trans == 'M') {
				for (iter = 0; iter < M; ++iter) {
					alphaval = *alpha * a[iter * M + iter];
					dscal(&SCALN, &alphaval, &c[iter], &incx);
				}
			} else if (*trans == 'V') {
				for (iter = 0; iter < M; ++iter) {
					alphaval = *alpha * a[iter];
					dscal(&SCALN, &alphaval, &c[iter], &incx);
				}
			}
		} else if (*side == 'R') {
			SCALN = M;
			incx = 1;
			if (*trans == 'M') {
				for (iter = 0; iter < N; ++iter) {
					alphaval = *alpha * a[iter * M + iter];
					dscal(&SCALN, &alphaval, &c[iter * M], &incx);
				}
			} else if (*trans == 'V') {
				for (iter = 0; iter < N; ++iter) {
					alphaval = *alpha * a[iter];
					dscal(&SCALN, &alphaval, &c[iter * M], &incx);
				}
			}
		}
	} else {
		BlasInt AXPYN;
		memset((void*) c, 0, M * N * sizeof(double));
		if (*side == 'L') {
			AXPYN = N;
			incx = M;
			if (*trans == 'M') {
				for (iter = 0; iter < M; ++iter) {
					alphaval = *alpha * a[iter * M + iter];
					daxpy(&AXPYN, &alphaval, &b[iter], &incx, &c[iter], &incx);
				}
			} else if (*trans == 'V') {
				for (iter = 0; iter < M; ++iter) {
					alphaval = *alpha * a[iter];
					daxpy(&AXPYN, &alphaval, &b[iter], &incx, &c[iter], &incx);
				}
			}
		} else if (*side == 'R') {
			AXPYN = M;
			incx = 1;
			if (*trans == 'M') {
				for (iter = 0; iter < N; ++iter) {
					alphaval = *alpha * a[iter * M + iter];
					daxpy(&AXPYN, &alphaval, &b[iter * M], &incx, &c[iter * M], &incx);
				}
			} else if (*trans == 'V') {
				for (iter = 0; iter < N; ++iter) {
					alphaval = *alpha * a[iter];
					daxpy(&AXPYN, &alphaval, &b[iter * M], &incx, &c[iter * M], &incx);
				}
			}
		}
	}
}

void ddiag(char* trans, BlasInt* m, BlasInt* n, double* a, double* b) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt iter;
	BlasInt MINMN = (M < N) ? M : N;
	if (*trans == 'V') {
		for (iter = 0; iter < MINMN; ++iter) {
			b[iter] = a[iter * M + iter];
		}
	} else if (*trans == 'M') {
		memset((void*) a, 0, M * N * sizeof(double));
		for (iter = 0; iter < MINMN; ++iter) {
			a[iter * M + iter] = b[iter];
		}
	}
}

double dtrac(BlasInt* m, BlasInt* n, double* a) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt iter;
	BlasInt MINMN = (M < N) ? M : N;
	double traceValue = 0;
	for (iter = 0; iter < MINMN; ++iter) {
		traceValue += a[iter * M + iter];
	}

	return traceValue;
}

double dvsum(BlasInt* n, double* x, BlasInt* incx) {

	double sum = 0;
	BlasInt iter;
	BlasInt N = *n;
	BlasInt INCX = *incx;

	for (iter = 0; iter < N; ++iter) {
		sum += x[iter * INCX];
	}

	return sum;
}

double dvprd(BlasInt* n, double* x, BlasInt* incx) {

	double product = 1;
	BlasInt iter;
	BlasInt N = *n;
	BlasInt INCX = *incx;

	for (iter = 0; iter < N; ++iter) {
		product *= x[iter * INCX];
	}

	return product;
}

void dgeva(char* trans, BlasInt* n, BlasInt* m, double* alpha, double* a, double* b) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt iter;
	BlasInt incx = 1;
	BlasInt incy;
	BlasInt AXPYN;

	if (*trans == 'C') {
		AXPYN = M;
		incy = 1;
		for (iter = 0; iter < N; ++iter) {
			daxpy(&AXPYN, alpha, a, &incx, &b[iter * M], &incy);
		}
	} else if (*trans == 'R') {
		AXPYN = N;
		incy = M;
		for (iter = 0; iter < M; ++iter) {
			daxpy(&AXPYN, alpha, a, &incx, &b[iter], &incy);
		}
	}
}

void dgevm(char* trans, BlasInt* m, BlasInt* n, double* alpha, double* a, double* b) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt incx;
	BlasInt iter;
	BlasInt SCALN;
	double alphaval;
	if (*trans == 'R') {
		SCALN = N;
		incx = M;
		for (iter = 0; iter < M; ++iter) {
			alphaval = *alpha * a[iter];
			dscal(&SCALN, &alphaval, &b[iter], &incx);
		}
	} else if (*trans == 'C') {
		SCALN = M;
		incx = 1;
		for (iter = 0; iter < N; ++iter) {
			alphaval = *alpha * a[iter];
			dscal(&SCALN, &alphaval, &b[iter * M], &incx);
		}
	}
}

void dgesu(char* trans, BlasInt* n, BlasInt* m, double* alpha, double* a, double* beta, double* b) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt iter;
	BlasInt incx;
	BlasInt incy = 1;
	BlasInt AXPYN;
	BlasInt SCALN;

	if (*trans == 'C') {
		AXPYN = M;
		SCALN = M;
		incx = 1;
		if (*beta == 0) {
			memset((void*) b, 0, SCALN * sizeof(double));
		} else {
			dscal(&SCALN, beta, b, &incy);
		}
		for (iter = 0; iter < N; ++iter) {
			daxpy(&AXPYN, alpha, &a[iter * M], &incx, b, &incy);
		}
	} else if (*trans == 'R') {
		AXPYN = N;
		SCALN = N;
		incx = M;
		if (*beta == 0) {
			memset((void* ) b, 0, SCALN * sizeof(double));
		} else {
			dscal(&SCALN, beta, b, &incy);
		}
		for (iter = 0; iter < M; ++iter) {
			daxpy(&AXPYN, alpha, &a[iter], &incx, b, &incy);
		}
	}
}
