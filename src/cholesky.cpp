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
	const BlasInt ione = 1;
	const BlasInt step = (*uplo == 'L') ? 1 : N;

	double cosine;
	double sine;
	BlasInt size = 0;
	double* data = NULL;

	*info = 0;
	for (BlasInt iter = 0; iter <= N - 1; ++iter) {
		data = &a[iter * N + iter];
		if ((*data == 0) && (work[iter] == 0)) {
			*info = 1;
			break;
		}
		drotg(data, &work[iter], &cosine, &sine);
		/*
		 * TODO: Change these to EPS equalities.
		 */
		if (*data < 0) {
			(*data) = - (*data);
			cosine = - cosine;
			sine = - sine;
		} else if (*data == 0) {
			*info = 1;
			break;
		}
		size = N - 1 - iter;
		drot(&size, &data[step], &step, &work[iter + 1], &ione, &cosine, &sine);
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
	const BlasInt step = (*uplo == 'L') ? 1 : N;

#if defined(USE_CHOLESKY_LINPACK)
	const char trans = (*uplo == 'L') ? 'N' : 'T';
	const char diag = 'N';
	dtrsv(uplo, &trans, &diag, n, a, n, work, &ione);

	*info = 0;
	double qs = dnrm2(n, work, &ione);
	qs = 1.0 - qs * qs;
	if (qs <= 0) {
		*info = 1;
	} else {
		BlasInt size = 0;
		double cosine;
		double sine;
		double* data = NULL;

		qs = sqrt(qs);
		memset((void*) work, 0, N * sizeof(work));
		for (BlasInt iter = N - 1; iter >= 0; --iter) {
			drotg(&qs, &work[iter], &cosine, &sine);
			if (qs <= 0) {
				qs = -qs;
				cosine = -cosine;
				sine = -sine;
			}
			data = &a[iter * (N + 1)];
			if (*data <= 0) {
				*info = 1;
				break;
			}
			size = N - iter;
			drot(&size, &work[iter], &ione, data, &step, &cosine, &sine);
			if (*data < 0) {
				qs = -1.0;
				dscal(&size, &qs, data, &step);
			/*
			 * TODO: Change these to EPS equalities.
			 */
			} else if (*data == 0) {
				*info = 1;
				break;
			}
		}
	}
#elif defined(USE_CHOLESKY_SIMD)
	double alpha = 1.0;
	double alphaPrev = alpha;
	double beta = 1.0;
	double betaPrev = beta;
	double ratio = 0;
	double invProduct = 0;
	double params[5];
	params[0] = -1.0;
	params[1] = 1.0;

	*info = 0;
	if (*uplo == 'U') {
		for (BlasInt iterI = 1; iterI <= N; ++iterI) {
			alphaPrev = alpha;
			betaPrev = beta;
			params[2] = - work[iterI - 1] / a[(iterI - 1) * (N + 1)];
			alpha = alphaPrev - params[2] * params[2];
			/*
			 * TODO: Change this to eps inequality.
			 */
			if (alpha <= 0) {
				*info = 1;
				break;
			}
			beta = sqrt(alpha);
			ratio = beta / betaPrev;
			invProduct = 1 / beta / betaPrev;
			params[3] = params[2] * invProduct;
			params[4] = alphaPrev * invProduct;
			a[(iterI - 1) * (N + 1)] *= ratio;
			BlasInt K = N - iterI;
			drotm(&K, &work[iterI], &ione, &a[(iterI - 1) + iterI * N ], &step, params);
		}
	} else {
		for (BlasInt iterI = 1; iterI <= N; ++iterI) {
			alphaPrev = alpha;
			betaPrev = beta;
			params[2] = - work[iterI - 1] / a[(iterI - 1) * (N + 1)];
			alpha = alphaPrev - params[2] * params[2];
			/*
			 * TODO: Change this to eps inequality.
			 */
			if (alpha <= 0) {
				*info = 1;
				break;
			}
			beta = sqrt(alpha);
			ratio = beta / betaPrev;
			invProduct = 1 / beta / betaPrev;
			params[3] = params[2] * invProduct;
			params[4] = alphaPrev * invProduct;
			a[(iterI - 1) * (N + 1)] *= ratio;
			BlasInt K = N - iterI;
			drotm(&K, &work[iterI], &ione, &a[iterI + (iterI - 1) * N ], &step, params);
		}
	}
#endif
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

void dchrk_sub(char* uplo, BlasInt* n, BlasInt* k, double* alpha, double* c,
			double* work, BlasInt* info) {
	for (BlasInt iter = 0; iter < *k; ++iter) {
		dchr_sub(uplo, n, alpha, c, &work[iter * *n], info);
		if (*info != 0) {
			break;
		}
	}
}

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

/*
 * TODO: Replace loops with BLAS calls, especially Givens rotations and
 * application.
 * TODO: Eliminate use of cvec and svec. If not possible, use lartv inside
 * double loop.
 */
void dchex(double* a, BlasInt* n, BlasInt* k, BlasInt* l, double* work,
		BlasInt* job) {
//	dchex_(a, &n, &n, &k, &l, xmat, &n, &nz, cvec, svec, &job);

	const BlasInt N = *n;
	const BlasInt L = *l;
	const BlasInt K = *k;
	const BlasInt lmk = L - K;
	const BlasInt lm1 = L - 1;
	const BlasInt ione = 1;
	const BlasInt imone = -1;
	double* cvec = work;
	double* svec = &work[N];

	if (*job == 1) {
		dcopy(&L, &a[lm1 * N], &imone, svec, &ione);

		for (BlasInt iterJ = lm1; iterJ >= K; --iterJ) {
			for (BlasInt iterI = 1; iterI <= iterJ; ++iterI) {
				a[iterJ * N + iterI - 1] = a[(iterJ - 1) * N + iterI - 1];
			}
			a[iterJ * N + iterJ] = 0;
		}

		BlasInt size = K - 1;
		dcopy(&size, &svec[L - K + 1], &imone, &a[(K - 1) * N], &ione);

		double t = svec[0];
		for (BlasInt iterI = 1; iterI <= lmk; ++iterI) {
			drotg(&svec[iterI], &t, &cvec[iterI - 1], &svec[iterI - 1]);
			t = svec[iterI];
		}

		a[(K - 1) * (N + 1)] = t;
		for (BlasInt iterJ = K + 1; iterJ <= N; ++iterJ) {
			BlasInt iterIL = (L - iterJ + 1 > 1)?(L - iterJ + 1):(1);
			for (BlasInt iterII = iterIL; iterII <= lmk; ++iterII) {
				BlasInt iterI = L - iterII;
				t = cvec[iterII - 1] * a[(iterJ - 1) * N + iterI - 1]
				  + svec[iterII - 1] * a[(iterJ - 1) * N + iterI];
				a[(iterJ - 1) * N + iterI] = cvec[iterII - 1]
				              * a[(iterJ - 1) * N + iterI]
				              - svec[iterII - 1]
				              * a[(iterJ - 1) * N + iterI - 1];
				a[(iterJ - 1) * N + iterI - 1] = t;
			}
		}
	} else if (*job == 2) {
		for (BlasInt iterI = 1; iterI <= K; ++iterI) {
			BlasInt iterII = lmk + iterI;
			svec[iterII - 1] = a[(K - 1) * N + iterI - 1];
		}

		for (BlasInt iterJ = K; iterJ <= lm1; ++iterJ) {
			for (BlasInt iterI = 1; iterI <= iterJ; ++iterI) {
				a[(iterJ - 1) * N + iterI - 1] = a[iterJ * N + iterI - 1];
			}
			BlasInt iterJJ = iterJ - K + 1;
			svec[iterJJ - 1] = a[iterJ * N + iterJ];
		}

		for (BlasInt iterI = 1; iterI <= K; ++iterI) {
			BlasInt iterII = lmk + iterI;
			a[(L - 1) * N + iterI - 1] = svec[iterII - 1];
		}

		memset((void *) &a[(L - 1) * N + K], 0, lmk * sizeof(double));

		for (BlasInt iterJ = K; iterJ <= N; ++iterJ) {
			if (iterJ != K) {
				BlasInt iterIU = (iterJ - 1 < L - 1)?(iterJ - 1):(L - 1);

				for (BlasInt iterI = K; iterI <= iterIU; ++iterI) {
					BlasInt iterII = iterI - K + 1;
					double t = cvec[iterII - 1] * a[(iterJ - 1) * N + iterI - 1]
					         + svec[iterII - 1] * a[(iterJ - 1) * N + iterI];
					a[(iterJ - 1) * N + iterI] = cvec[iterII - 1]
					                           * a[(iterJ - 1) * N + iterI]
					                           - svec[iterII - 1]
					                           * a[(iterJ - 1) * N + iterI - 1];
					a[(iterJ - 1) * N + iterI - 1] = t;
				}
			}

			if (iterJ < L) {
				BlasInt iterJJ = iterJ - K + 1;
				double t = svec[iterJJ - 1];
				drotg(&a[(iterJ - 1) * N + iterJ - 1], &t, &cvec[iterJJ - 1],
					&svec[iterJJ - 1]);
			}
		}
	}

	for (BlasInt iterJ = K; iterJ <= L; ++iterJ) {
		BlasInt iterJJ = N - iterJ + 1;
		double t = -1.0;
		dscal(&iterJJ, &t, &a[(iterJ - 1) * (N + 1)], &N);
	}
}

void dchd(double* a, BlasInt* n, BlasInt* k, double* work) {
	BlasInt job = 2;
	dchex(a, n, k, n, work, &job);
}

/*
 * TODO: This does not work because of continuous storage of a.
 */
//void dchdk(double* a, BlasInt* n, BlasInt* k, BlasInt* m, double* work) {
//	BlasInt job = 2;
//	for (BlasInt iter = 0; iter < *m; ++iter) {
//		BlasInt N = n - iter;
//		dchex(a, &N, &k[iter], &N, work, &job);
//	}
//}
