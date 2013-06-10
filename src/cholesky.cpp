/*
 * cholesky.h
 *
 *  Created on: Mar 15, 2012
 *      Author: igkiou
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cholesky.h"

void dchud_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info) {
	const BlasInt N = *n;
	double* wvec = work;
	double* cvec = &work[N];
	double* svec = &work[2 * N];
	const BlasInt ione = 1;
	const BlasInt step = (*uplo == 'L') ? 1 : N;
	BlasInt size = 0;
	double* tbuff = NULL;
	*info = 0;
	for (BlasInt iter = 0; iter < N - 1; ++iter) {
		tbuff = &a[iter * N + iter];
		if ((*tbuff == 0) && (wvec[iter] == 0)) {
			*info = 1;
			break;
		}
		drotg(tbuff, &wvec[iter], &cvec[iter], &svec[iter]);
		/*
		 * TODO: Change these to EPS equalities.
		 */
		if (*tbuff < 0) {
			(*tbuff) = - (*tbuff);
			cvec[iter] = - cvec[iter];
			svec[iter] = - svec[iter];
		} else if (*tbuff == 0) {
			*info = 1;
			break;
		}
		size = N - 1 - iter;
		if (size > 0) {
			drot(&size, &tbuff[step], &step, &wvec[iter + 1], &ione,
				&cvec[iter], &svec[iter]);
		}
	}

	tbuff = &a[(N - 1) * (N + 1)];
	if ((*info == 0) && ((*tbuff != 0) || (wvec[N - 1] != 0.0))) {
		drotg(tbuff, &wvec[N - 1], &cvec[N - 1], &svec[N - 1]);
		if (*tbuff < 0) {
			(*tbuff) = - (*tbuff);
			cvec[N - 1] = - cvec[N - 1];
			svec[N - 1] = - svec[N - 1];
		} else if (*tbuff == 0.0) {
			*info = 1;
		}
	} else {
		*info = 1;
	}
}

void dchud(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info) {
	const BlasInt ione = 1;
	dcopy(n, x, incx, work, &ione);
	dchud_sub(uplo, n, a, work, info);
}

void dchdd_sub(char* uplo, BlasInt* n, double* a, double* work, BlasInt* info) {
	const BlasInt N = *n;
	const BlasInt ione = 1;

	double* wvec = work;
	double* cvec = &work[N];
	double* svec = &work[2 * N];

	const char trans = (*uplo == 'L') ? 'N' : 'T';
	const char diag = 'N';
	dtrsv(uplo, &trans, &diag, n, a, n, wvec, &ione);

	*info = 0;
	double qs = dnrm2(n, wvec, &ione);
	qs = 1.0 - qs * qs;
	if (qs <= 0.0) {
		*info = 1;
	} else {
		qs = sqrt(qs);
		for (BlasInt iter = N - 1; iter >= 0; --iter) {
			drotg(&qs, &wvec[iter], &cvec[iter], &svec[iter]);
			if (qs <= 0.0) {
				qs = -qs;
				cvec[iter] = -cvec[iter];
				svec[iter] = -svec[iter];
			}
		}
		const BlasInt step = (*uplo == 'L') ? 1 : N;
		BlasInt size = 0;
		double* tbuff = NULL;
		memset((void*) wvec, 0, N * sizeof(wvec));
		for (BlasInt iter = N - 1; iter >= 0; --iter) {
			tbuff = &a[iter * (N + 1)];
			/*
			 * TODO: Change these to EPS equalities.
			 */
			if (*tbuff <= 0.0) {
				*info = 1;
				break;
			}
			size = N - iter;
			drot(&size, &wvec[iter], &ione, tbuff, &step, &cvec[iter],
				&svec[iter]);
			if (*tbuff < 0.0) {
				qs = -1.0;
				dscal(&size, &qs, tbuff, &step);
			} else if (*tbuff == 0.0) {
				*info = 1;
				break;
			}
		}
	}
}

void dchdd(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info) {
	const BlasInt ione = 1;
	dcopy(n, x, incx, work, &ione);
	dchdd_sub(uplo, n, a, work, info);
}

void dchr_sub(char* uplo, BlasInt* n, double* alpha, double* a, double* work,
			BlasInt* info) {
	if (*alpha > 0) {
		const BlasInt ione = 1;
		const double mult = sqrt(*alpha);
		dscal(n, &mult, work, &ione);
		dchud_sub(uplo, n, a, work, info);
	} else if (*alpha < 0) {
		const BlasInt ione = 1;
		const double mult = sqrt(- *alpha);
		dscal(n, &mult, work, &ione);
		dchdd_sub(uplo, n, a, work, info);
	} else {
		*info = 0;
	}
}

void dchr(char* uplo, BlasInt* n, double *alpha, double* x, BlasInt* incx,
		double* a, double* work, BlasInt* info) {
	const BlasInt ione = 1;
	dcopy(n, x, incx, work, &ione);
	dchr_sub(uplo, n, alpha, a, work, info);
}

void dchmv(char* uplo, BlasInt* n, double* a, double* x, BlasInt* incx) {
	const char diag = 'N';
	char transFirst;
	char transSecond;
	if (*uplo == 'U') {
		/* a is U, upper triangular, A * x = U' * (U * x) */
		transFirst = 'N';
		transSecond = 'T';
	} else if (*uplo == 'L') {
		/* a is L, lower triangular, A * x = L * (L' * x) */
		transFirst = 'T';
		transSecond = 'N';
	}
	dtrmv(uplo, &transFirst, &diag, n, a, n, x, incx);
	dtrmv(uplo, &transSecond, &diag, n, a, n, x, incx);
}

/*
 * TODO: May make sense to provide option(s) for
 * 1) work being the size of a plus 2; or
 * 2) custom dchud and dchdd working directly on v and work of size 2*n; or
 * 3) only avoid the multiple rescalings in dchr, and otherwise keep the same.
 */
void dchrk(char* uplo, char* trans, BlasInt* n, BlasInt* k, double* alpha,
		double* a, double* c, double* work, BlasInt* info) {
	if (*trans == 'N') {
		BlasInt incx = 1;
		for (BlasInt iter = 0; iter < *k; ++iter) {
			dchr(uplo, n, alpha, &a[iter * *n], &incx, c, work, info);
			if (*info != 0) {
				break;
			}
		}
	} else if (*trans == 'L') {
		BlasInt incx = *k;
		for (BlasInt iter = 0; iter < *k; ++iter) {
			dchr(uplo, n, alpha, &a[iter], &incx, c, work, info);
			if (*info != 0) {
				break;
			}
		}
	}
}

void dchmm(char* side, char* uplo, BlasInt* m, BlasInt* n, double* alpha,
		double* a, double* b) {
	const char diag = 'N';
	const double done = 1.0;
	char transFirst;
	char transSecond;
	BlasInt* lda;
	if (*side == 'L') {
		lda = m;
		if (*uplo == 'U') {
			/* a is U, upper triangular, A * B = U' * (U * B) */
			transFirst = 'N';
			transSecond = 'T';
		} else if (*uplo == 'L') {
			/* a is L, lower triangular, A * B = L * (L' * B) */
			transFirst = 'T';
			transSecond = 'N';
		}
	} else if (*side == 'R') {
		lda = n;
		if (*uplo == 'U') {
			/* a is U, upper triangular, B * A = (B * U') * U */
			transFirst = 'T';
			transSecond = 'N';
		} else if (*uplo == 'L') {
			/* a is L, lower triangular, B * A = (B * L) * L' */
			transFirst = 'N';
			transSecond = 'T';
		}
	}
	dtrmm(side, uplo, &transFirst, &diag, m, n, &done, a, lda, b, m);
	dtrmm(side, uplo, &transSecond, &diag, m, n, alpha, a, lda, b, m);
}
