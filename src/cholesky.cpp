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
	BlasInt size = N;
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
		if (--size > 0) {
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

	const char trans[2] = {(*uplo == 'L') ? 'N' : 'T', 0};
	const char diag = 'N';
	dtrsv(uplo, trans, &diag, n, a, n, wvec, &ione);

	*info = 0;
	double qs = dnrm2(n, wvec, &ione);
	qs = 0.0, 1.0 - qs * qs;
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
		BlasInt size = N;
		BlasInt npos;
		BlasInt* flind = NULL;
		double* tbuff = NULL;
		memset((void*) wvec, 0, N * sizeof(wvec));
		for (BlasInt iter = N - 1; iter >= 0; --iter) {
			tbuff = &a[iter * (N + 1)];
			--size;
			if (*tbuff <= 0.0) {
				*info = 1;
				break;
			}
			drot(&size, &wvec[iter], &ione, tbuff, &step, &cvec[iter], &svec[iter]);
			if (*tbuff < 0.0) {
				if (flind == NULL) {
					flind = (BlasInt*) malloc(N * sizeof(BlasInt));
					npos = N;
				}
				flind[--npos] = iter;
				qs = -1.0;
				dscal(&size, &qs, tbuff, &step);
			} else if (*tbuff == 0.0) {
				*info = 1;
				break;
			}
		}
		if (flind != NULL) {
			free((void*) flind);
		}
	}
}

void dchdd(char* uplo, BlasInt* n, double* x, BlasInt* incx, double* a,
		double* work, BlasInt* info) {
	const BlasInt ione = 1;
	dcopy(n, x, incx, work, &ione);
	dchdd_sub(uplo, n, a, work, info);
}

//void dchdd(char* uplo, BlasInt* pN, double* a, double* work, BlasInt* info) {
//
//	BlasInt N = *pN;
//	double* cvec = work;
//	double* svec = &work[N];
//	double* wvec = &work[2 * N];
//	BlasInt incx = 1;
//	BlasInt retcode = 0;
//	char trans = (*uplo == 'L') ? 'N' : 'T';
//	char diag = 'N';
//	dtrsv(uplo, &trans, &diag, &N, a, &N, wvec, &incx);
//	double temp = dnrm2(&N, wvec, &incx);
//	double qs = (1.0 - temp) * (1.0 + temp);
//	int* flind = NULL;
//	BlasInt npos;
//	if (qs <= 0) {
//		retcode = 1;
//	} else {
//		qs = sqrt(qs);
//		for (BlasInt iter = N - 1; iter >= 0; --iter) {
//			drotg(&qs, &wvec[iter], &cvec[iter], &svec[iter]);
//			if (qs < 0) {
//				qs = - qs;
//				cvec[iter] = - cvec[iter];
//				svec[iter] = - svec[iter];
//			}
//		}


//	    fillVec(wkvec,n,0.0);
//	    stp=islower?1:lmat.stride;
//	    for (i=n-1,sz=0,tbuff=lmat.buff+((n-1)*(lmat.stride+1)); i>=0; i--) {
//	      /* BAD: Slower for upper triangular! */
//	      sz++;
//	      if (*tbuff<=0.0) {
//		retcode=1; break;
//	      }
//	      BLASFUNC(drot) (&sz,wkvec+i,&ione,tbuff,&stp,cvec+i,svec+i);
//	      /* Do not want negative elements on diagonal */
//	      if (*tbuff<0.0) {
//		if (flind==0) {
//		  /* Does this ever happen?
//		     Allocate 'flind'. Size n, to make sure. Grows from the right */
//		  flind=(int*) mxMalloc(n*sizeof(int));
//		  npos=n;
//		}
//		flind[--npos]=i;
//		qs=-1.0;
//		BLASFUNC(dscal) (&sz,&qs,tbuff,&stp);
//	      } else if (*tbuff==0.0) {
//		retcode=1; break;
//	      }
//	      tbuff-=(lmat.stride+1);
//	    }
	    /* NOTE: Should have v in 'wkvec' now */

//		if (*uplo == 'U') {
//
//		}
//		memset(wvec, 0, N * sizeof(wvec));
//		BlasInt step = (*uplo == 'L') ? 1 : N;
//		double *tbuff = &a[(N - 1) * (N + 1)];
//		BlasInt size = 0;
//		for (BlasInt iter = N - 1; iter >= 0; --iter) {
//			tbuff = &a[iter * (N + 1)];
//			if (*tbuff <= 0) {
//				retcode = 1;
//				break;
//			}
//			++size;
//			drot(&size, &wvec[iter], &incx, tbuff, &step, &cvec[iter], &svec[iter]);
//			if (*tbuff < 0) {
//				if (flind == NULL) {
//					flind = (int *) malloc(N * sizeof(int));
//					npos = N;
//				}
//				flind[--npos] = iter;
//				qs = - 1.0;
//				dscal(&size, &qs, tbuff, &step);
//			} else
//			if (*tbuff == 0) {
//				retcode = 1;
//				break;
//			}
//			tbuff -= (N + 1);
//		}
//	}
//	if (flind != NULL) {
//		free(flind);
//	}
//	*info = retcode;
//}

//void dchr(char* uplo, BlasInt* N, double* alpha, double* x, BlasInt* incx, \
//		double* a, BlasInt* lda, double* work, BlasInt* lwork, BlasInt* info) {
//
//	if (*lda != *N) {
//		*info = - 7;
//		return;
//	}
//	if (*lwork == -1) {
//		*work = (double) 3 * (*N);
//		*info = 0;
//		return;
//	} else if (*lwork < 3 * *N) {
//		*info = - 9;
//		return;
//	}
//
//	double* wvec = &work[2 * *N];
//	double mult = sqrt(*alpha);
//	BlasInt incy = 1;
//	dcopy(N, x, incx, wvec, &incy);
//	dscal(N, &mult, wvec, &incy);
//	if (*alpha > 0) {
//		dchud(uplo, N, a, work, info);
//	} else if (*alpha < 0){
//		dchdd(uplo, N, a, work, info);
//	} else {
//		*info = 0;
//	}
//}
//
//void dchmv(char* uplo, BlasInt* N, double* a, BlasInt* lda, double* x, \
//			BlasInt* incx) {
//	char diag = 'N';
//	if (*uplo == 'U') {
//		/* a is U, upper triangular, A * x = U' * (U * x) */
//		char trans = 'N';
//		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
//		trans = 'T';
//		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
//	} else if (*uplo == 'L') {
//		/* a is L, lower triangular, A * x = L * (L' * x) */
//		char trans = 'T';
//		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
//		trans = 'N';
//		dtrmv(uplo, &trans, &diag, N, a, lda, x, incx);
//	}
//}
//
//void dchrk(char* uplo, char* trans, BlasInt* N, BlasInt* K, double* alpha, \
//			double* a, BlasInt* lda, double* c, BlasInt* ldc, double* work, \
//			BlasInt* lwork, BlasInt* info) {
//
//	if ((*trans == 'N') && (*lda != *N)) {
//		*info = - 7;
//		return;
//	} else if ((*trans == 'T') && (*lda != *K)) {
//		*info = - 7;
//		return;
//	}
//	if (*ldc != *N) {
//		*info = - 9;
//		return;
//	}
//	if (*lwork == -1) {
//		*work = (double) 3 * *N;
//		*info = 0;
//		return;
//	} else if (*lwork < 3 * *N) {
//		*info = - 11;
//		return;
//	}
//
//	if (*trans == 'N') {
//		BlasInt incx = 1;
//		for (BlasInt iter = 0; iter < *K; ++iter) {
//			dchr(uplo, N, alpha, &a[iter * *N], &incx, c, lda, work, lwork, info);
//		}
//	} else if (*trans == 'L') {
//		BlasInt incx = *K;
//		for (BlasInt iter = 0; iter < *K; ++iter) {
//			dchr(uplo, N, alpha, &a[iter], &incx, c, lda, work, lwork, info);
//		}
//	}
//}
//
//void dchmm(char* uplo, char* side, BlasInt* M, BlasInt* N, double* alpha, \
//			double* a, BlasInt* lda, double* b, BlasInt* ldb) {
//	char diag = 'N';
//	double done = 1.0;
//	if (*side == 'L') {
//		if (*uplo == 'U') {
//			/* a is U, upper triangular, A * B = U' * (U * B) */
//			char trans = 'N';
//			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
//			trans = 'T';
//			dtrmm(side, uplo, &trans, &diag, N, N, alpha, a, lda, b, ldb);
//		} else if (*uplo == 'L') {
//			/* a is L, lower triangular, A * B = L * (L' * B) */
//			char trans = 'T';
//			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
//			trans = 'N';
//			dtrmm(side, uplo, &trans, &diag, M, N, alpha, a, lda, b, ldb);
//		}
//	} else if (*side == 'R') {
//		if (*uplo == 'U') {
//			/* a is U, upper triangular, B * A = (B * U') * U */
//			char trans = 'T';
//			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
//			trans = 'N';
//			dtrmm(side, uplo, &trans, &diag, N, N, alpha, a, lda, b, ldb);
//		} else if (*uplo == 'L') {
//			/* a is L, lower triangular, B * A = (B * L) * L' */
//			char trans = 'N';
//			dtrmm(side, uplo, &trans, &diag, M, N, &done, a, lda, b, ldb);
//			trans = 'T';
//			dtrmm(side, uplo, &trans, &diag, M, N, alpha, a, lda, b, ldb);
//		}
//	}
//}
