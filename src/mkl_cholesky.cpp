/*
 * cholesky.h
 *
 *  Created on: Mar 15, 2012
 *      Author: igkiou
 */

#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mkl_cholesky.h"

void dchud(char *uplo, MKL_INT *pN, double *a, double *work, MKL_INT *info) {

	MKL_INT N = *pN;
	double *cvec = work;
	double *svec = &work[N];
	double *wvec = &work[2 * N];
	MKL_INT incx = 1;
	MKL_INT retcode = 0;
	MKL_INT step = (*uplo == 'L') ? 1 : N;
	double *tbuff;
	MKL_INT size;
	for (MKL_INT iter = 0; iter < N - 1; ++iter) {
		tbuff = &a[iter * N + iter];
		if ((*tbuff == 0) && (wvec[iter] == 0)) {
			retcode = 1;
			break;
		}
		drotg(tbuff, &wvec[iter], &cvec[iter], &svec[iter]);
		if (*tbuff < 0) {
			*tbuff= - *tbuff;
			cvec[iter] = - cvec[iter];
			svec[iter] = - svec[iter];
		} else if (*tbuff == 0) {
			retcode = 1;
			break;
		}
		size = N - iter;
		drot(&size, &tbuff[step], &step, &wvec[iter + 1], &incx, &cvec[iter], \
			&svec[iter]);
	}

	tbuff = &a[(N - 1) * (N + 1)];
	if ((retcode == 0) && ((*tbuff != 0) || (wvec[N - 1] != 0))) {
		drotg(tbuff, &wvec[N - 1], &cvec[N - 1], &svec[N - 1]);
		if (*tbuff < 0) {
			*tbuff = - *tbuff;
			cvec[N - 1] = - cvec[N - 1];
			svec[N - 1] = - svec[N - 1];
		} else if (*tbuff == 0) {
			retcode = 1;
		}
	} else {
		retcode = 1;
	}
	*info = retcode;
}

void dchdd(char *uplo, MKL_INT *pN, double *a, double *work, MKL_INT *info) {

	MKL_INT N = *pN;
	double *cvec = work;
	double *svec = &work[N];
	double *wvec = &work[2 * N];
	MKL_INT incx = 1;
	MKL_INT retcode = 0;
	char trans = (*uplo == 'L') ? 'N' : 'T';
	char diag = 'N';
	dtrsv(uplo, &trans, &diag, &N, a, &N, wvec, &incx);
	double temp = dnrm2(&N, wvec, &incx);
	double qs = 1.0 - temp * temp;
	int *flind;
	int npos;
	if (qs <= 0) {
		retcode = 1;
	} else {
		qs = sqrt(qs);
		for (MKL_INT iter = N - 1; iter >=0; --iter) {
			drotg(&qs, &wvec[iter], &cvec[iter], &svec[iter]);
			if (qs < 0) {
				qs = - qs;
				cvec[iter] = - cvec[iter];
				svec[iter] = - svec[iter];
			}
		}

		memset(wvec, 0, N * sizeof(wvec));
		MKL_INT step = (*uplo == 'L') ? 1 : N;
		double *tbuff;
		MKL_INT size;
		for (MKL_INT iter = N - 1; iter >= 0; --iter) {
			size = iter;
			tbuff = &a[iter * (N + 1)];
			if (*tbuff <= 0) {
				retcode = 1;
				break;
			}
			drot(&size, &wvec[iter], &incx, tbuff, &step, &cvec[iter], &svec[iter]);
			if (*tbuff <= 0) {
				if (flind == 0) {
					flind = (int *) malloc(N * sizeof(int));
					npos = N;
				}
				flind[--npos] = iter;
				qs = - 1.0;
				dscal(&size, &qs, tbuff, &step);
			} else if (*tbuff == 0) {
				retcode = 1;
				break;
			}
		}
	}
	if (flind != 0) {
		free(flind);
	}
	*info = retcode;
}

void dchr(char *uplo, MKL_INT *N, double *alpha, double *x, MKL_INT *incx, \
		double *a, MKL_INT *lda, double *work, MKL_INT *lwork, MKL_INT *info) {

	if (*lda != *N) {
		*info = - 7;
		return;
	}
	if (*lwork == -1) {
		*work = (double) 3 * (*N);
		*info = 0;
		return;
	} else if (*lwork < 3 * *N) {
		*info = - 9;
		return;
	}

	double *wvec = &work[2 * *N];
	double mult = sqrt(*alpha);
	MKL_INT incy = 1;
	dcopy(N, x, incx, wvec, &incy);
	dscal(N, &mult, wvec, &incy);
	if (*alpha > 0) {
		dchud(uplo, N, a, work, info);
	} else if (*alpha < 0){
		dchdd(uplo, N, a, work, info);
	} else {
		*info = 0;
	}
}

void dchmv(char *uplo, MKL_INT *N, double *a, MKL_INT *lda, double *x, \
			MKL_INT *incx) {
	char diag = 'N';
	if (*uplo == 'U') {
		/* a is U, upper triangular, A * x = U' * (U * x) */
		char trans = 'N';
		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
		trans = 'T';
		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
	} else if (*uplo == 'L') {
		/* a is L, lower triangular, A * x = L * (L' * x) */
		char trans = 'T';
		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
		trans = 'N';
		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
	}
}

void dchrk(char *uplo, char *trans, MKL_INT *N, MKL_INT *K, double *alpha, \
			double *a, MKL_INT *lda, double *c, MKL_INT *ldc, double *work, \
			MKL_INT *lwork, MKL_INT *info) {

	if ((*trans == 'N') && (*lda != *N)) {
		*info = - 7;
		return;
	} else if ((*trans == 'T') && (*lda != *K)) {
		*info = - 7;
		return;
	}
	if (*ldc != *N) {
		*info = - 9;
		return;
	}
	if (*lwork == -1) {
		*work = (double) 3 * *N;
		*info = 0;
		return;
	} else if (*lwork < 3 * *N) {
		*info = - 11;
		return;
	}

	if (*trans == 'N') {
		MKL_INT incx = 1;
		for (MKL_INT iter = 0; iter < *K; ++iter) {
			dchr(uplo, N, alpha, &a[iter * *N], &incx, c, lda, work, lwork, info);
		}
	} else if (*trans == 'L') {
		MKL_INT incx = *K;
		for (MKL_INT iter = 0; iter < *K; ++iter) {
			dchr(uplo, N, alpha, &a[iter], &incx, c, lda, work, lwork, info);
		}
	}
}

void dchmm(char *uplo, char *side, MKL_INT *M, MKL_INT *N, double *alpha, \
			double *a, MKL_INT *lda, double *b, MKL_INT *ldb) {
	char diag = 'N';
	double done = 1.0;
	if (*side == 'L') {
		if (*uplo == 'U') {
			/* a is U, upper triangular, A * B = U' * (U * B) */
			char trans = 'N';
			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
			trans = 'T';
			dtrmm(side, uplo, &trans, &diag, N, N, alpha, a, lda, b, ldb);
		} else if (*uplo == 'L') {
			/* a is L, lower triangular, A * B = L * (L' * B) */
			char trans = 'T';
			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
			trans = 'N';
			dtrmm(side, uplo, &trans, &diag, M, N, alpha, a, lda, b, ldb);
		}
	} else if (*side == 'R') {
		if (*uplo == 'U') {
			/* a is U, upper triangular, B * A = (B * U') * U */
			char trans = 'T';
			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
			trans = 'N';
			dtrmm(side, uplo, &trans, &diag, N, N, alpha, a, lda, b, ldb);
		} else if (*uplo == 'L') {
			/* a is L, lower triangular, B * A = (B * L) * L' */
			char trans = 'N';
			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
			trans = 'T';
			dtrmm(side, uplo, &trans, &diag, M, N, alpha, a, lda, b, ldb);
		}
	}
}
