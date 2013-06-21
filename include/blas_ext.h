/*
 * mkl_ext//mkl_ext/include/blas_helper.h/blas_helper.h
 *
 *  Created on: Jun 6, 2013
 *      Author: igkiou
 */

#ifndef BLAS_EXT_H_
#define BLAS_EXT_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "blas_header.h"
#include "rng.h"
#include "cholesky.h"
#include "qr.h"
#include "arpack.h"
#include "propack.h"
#include "utils.h"

/*
 * BLAS helper functions.
 */
inline void blasXerbla(const char* srname, const int* info, const int lsrname) {
	xerbla(srname, info, lsrname);
}
inline int blasLsame(const char* ca, const char* cb, const BlasInt lca,
			const BlasInt lcb) {
	return lsame(ca, cb, lca, lcb);
}

/*
 * BLAS Level 1.
 */
inline double blasDasum(const BlasInt* n, const double* x,
						const BlasInt* incx) {
	return dasum(n, x, incx);
}
inline void blasDaxpy(const BlasInt* n, const double* alpha, const double* x,
					const BlasInt* incx, double* y, const BlasInt* incy) {
	daxpy(n, alpha, x, incx, y, incy);
}
inline void blasDaxpby(const BlasInt* n, const double* alpha, const double* x,
					const BlasInt* incx, const double* beta, double* y,
					const BlasInt* incy) {
	daxpby(n, alpha, x, incx, beta, y, incy);
}
inline void blasDaxpyi(const BlasInt* nz, const double* a, const double* x,
					const BlasInt* indx, double* y) {
	daxpyi(nz, a, x, indx, y);
}
inline void blasDcopy(const BlasInt* n, const double* x, const BlasInt* incx,
					double* y, const BlasInt* incy) {
	dcopy(n, x, incx, y, incy);
}
inline double blasDdot(const BlasInt* n, const double* x, const BlasInt* incx,
					const double* y, const BlasInt* incy) {
	return ddot(n, x, incx, y, incy);
}
inline double blasDsdot(const BlasInt* n, const float* x, const BlasInt* incx,
						const float* y, const BlasInt* incy) {
	return dsdot(n, x, incx, y, incy);
}
inline double blasDdoti(const BlasInt* nz, const double* x, const BlasInt* indx,
					const double* y) {
	return ddoti(nz, x, indx, y);
}
inline void blasDgthr(const BlasInt* nz, const double* y, double* x,
					const BlasInt* indx) {
	dgthr(nz, y, x, indx);
}
inline void blasDgthrz(const BlasInt* nz, double* y, double* x,
					const BlasInt* indx) {
	dgthrz(nz, y, x, indx);
}
inline double blasDnrm2(const BlasInt* n, const double* x,
						const BlasInt* incx) {
	return dnrm2(n, x, incx);
}
inline void blasDrot(const BlasInt* n, double* x, const BlasInt* incx,
					double* y, const BlasInt* incy, const double* c,
					const double* s) {
	drot(n, x, incx, y, incy, c, s);
}
inline void blasDrotg(double* a, double* b, double* c, double* s) {
	drotg(a, b, c, s);
}
inline void blasDroti(const BlasInt* nz, double* x, const BlasInt* indx,
					double* y, const double* c, const double* s) {
	droti(nz, x, indx, y, c, s);
}
inline void blasDrotm(const BlasInt* n, double* x, const BlasInt* incx,
					double* y, const BlasInt* incy, const double* param) {
	drotm(n, x, incx, y, incy, param);
}
inline void blasDrotmg(double* d1, double* d2, double* x1, const double* y1,
					double* param) {
	drotmg(d1, d2, x1, y1, param);
}
inline void blasDscal(const BlasInt* n, const double* a, double* x,
					const BlasInt* incx) {
	dscal(n, a, x, incx);
}
inline void blasDsctr(const BlasInt* nz, const double* x, const BlasInt* indx,
					double* y) {
	dsctr(nz, x, indx, y);
}
inline void blasDswap(const BlasInt* n, double* x, const BlasInt* incx,
					double* y, const BlasInt* incy) {
	dswap(n, x, incx, y, incy);
}
inline BlasInt blasIdamax(const BlasInt* n, const double* x,
						const BlasInt* incx) {
	return idamax(n, x, incx);
}
inline BlasInt blasIdamin(const BlasInt* n, const double* x,
						const BlasInt* incx) {
	return idamin(n, x, incx);
}

/*
 * BLAS Level 2.
 */
inline void blasDgbmv(const char* trans, const BlasInt* m, const BlasInt* n,
					const BlasInt* kl, const BlasInt* ku, const double* alpha,
					const double* a, const BlasInt* lda, const double* x,
					const BlasInt* incx, const double* beta, double* y,
					const BlasInt* incy) {
	dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}
inline void blasDgemv(const char* trans, const BlasInt* m, const BlasInt* n,
					const double* alpha, const double* a, const BlasInt* lda,
					const double* x, const BlasInt* incx, const double* beta,
					double* y, const BlasInt* incy) {
	dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void blasDger(const BlasInt* m, const BlasInt* n, const double* alpha,
					const double* x, const BlasInt* incx, const double* y,
					const BlasInt* incy, double* a, const BlasInt* lda) {
	dger(m, n, alpha, x, incx, y, incy, a, lda);
}
inline void blasDsbmv(const char* uplo, const BlasInt* n, const BlasInt* k,
					const double* alpha, const double* a, const BlasInt* lda,
					const double* x, const BlasInt* incx, const double* beta,
					double* y, const BlasInt* incy) {
	dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
}
inline void blasDspmv(const char* uplo, const BlasInt* n, const double* alpha,
					const double* ap, const double* x, const BlasInt* incx,
					const double* beta, double* y, const BlasInt* incy) {
	dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy);
}
inline void blasDspr(const char* uplo, const BlasInt* n, const double* alpha,
					const double* x, const BlasInt* incx, double* ap) {
	dspr(uplo, n, alpha, x, incx, ap);
}
inline void blasDspr2(const char* uplo, const BlasInt* n, const double* alpha,
					const double* x, const BlasInt* incx, const double* y,
					const BlasInt* incy, double* ap) {
	dspr2(uplo, n, alpha, x, incx, y, incy, ap);
}
inline void blasDsymv(const char* uplo, const BlasInt* n, const double* alpha,
					const double* a, const BlasInt* lda, const double* x,
					const BlasInt* incx, const double* beta, double* y,
					const BlasInt* incy) {
	dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline void blasDsyr(const char* uplo, const BlasInt* n, const double* alpha,
					const double* x, const BlasInt* incx, double* a,
					const BlasInt* lda) {
	dsyr(uplo, n, alpha, x, incx, a, lda);
}
inline void blasDsyr2(const char* uplo, const BlasInt* n, const double* alpha,
					const double* x, const BlasInt* incx, const double* y,
					const BlasInt* incy, double* a, const BlasInt* lda) {
	dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda);
}
inline void blasDtbmv(const char* uplo, const char* trans, const char* diag,
					const BlasInt* n, const BlasInt* k, const double* a,
					const BlasInt* lda, double* x, const BlasInt* incx) {
	dtbmv(uplo, trans, diag, n, k, a, lda, x, incx);
}
inline void blasDtbsv(const char* uplo, const char* trans, const char* diag,
					const BlasInt* n, const BlasInt* k, const double* a,
					const BlasInt* lda, double* x, const BlasInt* incx) {
	dtbsv(uplo, trans, diag, n, k, a, lda, x, incx);
}
inline void blasDtpmv(const char* uplo, const char* trans, const char* diag,
					const BlasInt* n, const double* ap, double* x,
					const BlasInt* incx) {
	dtpmv(uplo, trans, diag, n, ap, x, incx);
}
inline void blasDtpsv(const char* uplo, const char* trans, const char* diag,
					const BlasInt* n, const double* ap, double* x,
					const BlasInt* incx) {
	dtpsv(uplo, trans, diag, n, ap, x, incx);
}
inline void blasDtrmv(const char* uplo, const char* transa, const char* diag,
					const BlasInt* n, const double* a, const BlasInt* lda,
					double* b, const BlasInt* incx) {
	dtrmv(uplo, transa, diag, n, a, lda, b, incx);
}
inline void blasDtrsv(const char* uplo, const char* trans, const char* diag,
					const BlasInt* n, const double* a, const BlasInt* lda,
					double* x, const BlasInt* incx) {
	dtrsv(uplo, trans, diag, n, a, lda, x, incx);
}
inline void blasDgem2vu(const BlasInt* m, const BlasInt* n, const double* alpha,
						const double* a, const BlasInt* lda, const double* x1,
						const BlasInt* incx1, const double* x2,
						const BlasInt* incx2, const double* beta, double* y1,
						const BlasInt* incy1, double* y2,
						const BlasInt* incy2) {
	dgem2vu(m, n, alpha, a, lda, x1, incx1, x2, incx2, beta, y1, incy1, y2,
			incy2);
}

/*
 * BLAS Level 3.
 */
inline void blasDgemm(const char* transa, const char* transb, const BlasInt* m,
					const BlasInt* n, const BlasInt* k, const double* alpha,
					const double* a, const BlasInt* lda, const double* b,
					const BlasInt* ldb, const double* beta, double* c,
					const BlasInt* ldc) {
	dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void blasDsymm(const char* side, const char* uplo, const BlasInt* m,
					const BlasInt* n, const double* alpha, const double* a,
					const BlasInt* lda, const double* b, const BlasInt* ldb,
					const double* beta, double* c, const BlasInt* ldc) {
	dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void blasDsyr2k(const char* uplo, const char* trans, const BlasInt* n,
					const BlasInt* k, const double* alpha, const double* a,
					const BlasInt* lda, const double* b, const BlasInt* ldb,
					const double* beta, double* c, const BlasInt* ldc) {
	dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void blasDsyrk(const char* uplo, const char* trans, const BlasInt* n,
					const BlasInt* k, const double* alpha, const double* a,
					const BlasInt* lda, const double* beta, double* c,
					const BlasInt* ldc) {
	dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
}
inline void blasDtrmm(const char* side, const char* uplo, const char* transa,
					const char* diag, const BlasInt* m, const BlasInt* n,
					const double* alpha, const double* a, const BlasInt* lda,
					double* b, const BlasInt* ldb) {
	dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}
inline void blasDtrsm(const char* side, const char* uplo, const char* transa,
					const char* diag, const BlasInt* m, const BlasInt* n,
					const double* alpha, const double* a, const BlasInt* lda,
					double* b, const BlasInt* ldb) {
	dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

/*
 * LAPACK selection.
 */
inline void blasDgelsy(BlasInt* m, BlasInt* n, BlasInt* nrhs, double* a,
					BlasInt* lda, double* b, BlasInt* ldb, BlasInt* jpvt,
					double* rcond, BlasInt* rank, double* work, BlasInt* lwork,
					BlasInt* info) {
	dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
inline void blasDgesdd(char* jobz, BlasInt* m, BlasInt* n, double* a,
					BlasInt* lda, double* s, double* u, BlasInt* ldu,
					double* vt, BlasInt* ldvt, double* work, BlasInt* lwork,
					BlasInt* iwork, BlasInt* info) {
	dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline void blasDgesvd(char* jobu, char* jobvt, BlasInt* m, BlasInt* n,
					double* a, BlasInt* lda, double* s, double* u, BlasInt* ldu,
					double* vt, BlasInt* ldvt, double* work, BlasInt* lwork,
					BlasInt* info) {
	dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
inline void blasDlacpy(char* uplo, BlasInt* m, BlasInt* n, double* a,
					BlasInt* lda, double* b, BlasInt* ldb) {
	dlacpy(uplo, m, n, a, lda, b, ldb);
}
inline double blasDlamch(char* cmach) {
	return dlamch(cmach);
}
inline double blasDlange(char* norm, BlasInt* m, BlasInt* n, double* a,
						BlasInt* lda, double* work) {
	return dlange(norm, m, n, a, lda, work);
}
inline double blasDlansy(char* norm, char* uplo, BlasInt* n, double* a,
						BlasInt* lda, double* work) {
	return dlansy(norm, uplo, n, a, lda, work);
}
inline void blasDpotrf(char* uplo, BlasInt* n, double* a, BlasInt* lda,
					BlasInt* info) {
	dpotrf(uplo, n, a, lda, info);
}
inline void blasDpotrs(char* uplo, BlasInt* n, BlasInt* nrhs, double* a,
					BlasInt* lda, double* b, BlasInt* ldb, BlasInt* info) {
	dpotrs(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline void blasDsyev(char* jobz, char* uplo, BlasInt* n, double* a,
					BlasInt* lda, double* w, double* work, BlasInt* lwork,
					BlasInt* info) {
	dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
}
inline void blasDsyevr(char* jobz, char* range, char* uplo, BlasInt* n,
					double* a, BlasInt* lda, double* vl, double* vu,
					BlasInt* il, BlasInt* iu, double* abstol, BlasInt* m,
					double* w, double* z, BlasInt* ldz, BlasInt* isuppz,
					double* work, BlasInt* lwork, BlasInt* iwork,
					BlasInt* liwork, BlasInt* info) {
	dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
		isuppz, work, lwork, iwork, liwork, info);
}

/*
 * Cholesky update.
 */
inline void blasDchud(char* uplo, BlasInt* n, double* x, BlasInt* incx,
					double* a, double* work, BlasInt* info) {
	dchud(uplo, n, x, incx, a, work, info);
}

inline void blasDchdd(char* uplo, BlasInt* n, double* x, BlasInt* incx,
					double* a, double* work, BlasInt* info) {
	dchdd(uplo, n, x, incx, a, work, info);
}

inline void blasDchr(char* uplo, BlasInt* n, double* alpha, double* x,
					BlasInt* incx, double* a, double* work, BlasInt* info) {
	dchr(uplo, n, alpha, x, incx, a, work, info);
}

inline void blasDchmv(char* uplo, BlasInt* n, double* a, double* x,
					BlasInt* incx) {
	dchmv(uplo, n, a, x, incx);
}

inline void blasDchrk(char* uplo, char* trans, BlasInt* n, BlasInt* k,
					double* alpha, double* a, double* c, double* work,
					BlasInt* info) {
	dchrk(uplo, trans, n, k, alpha, a, c, work, info);
}

inline void blasDchmm(char* side, char* uplo, BlasInt* m, BlasInt* n,
					double* alpha, double* a, double* b) {
	dchmm(side, uplo, m, n, alpha, a, b);
}

inline void blasDchd(double* a, BlasInt* n, BlasInt* k, double* work) {
	dchd(a, n, k, work);
}

/*
 * QR update.
 */

/*
 * ARPACK.
 */

/*
 * PROPACK.
 */

/*
 * Utilities.
 */
inline void blasDdimm(char* side, char* trans, BlasInt* m, BlasInt* n,
					double* alpha, double* a, double* b, double* beta,
					double* c) {
	ddimm(side, trans, m, n, alpha, a, b, beta, c);
}
inline void blasDdiag(char* trans, BlasInt* m, BlasInt* n, double* a,
					double* b) {
	ddiag(trans, m, n, a, b);
}
inline double blasDtrac(BlasInt* m, BlasInt* n, double* a) {
	return dtrac(m, n, a);
}
inline double blasDvsum(BlasInt* n, double* x, BlasInt* incx) {
	return dvsum(n, x, incx);
}
inline double blasDvprd(BlasInt* n, double* x, BlasInt* incx) {
	return dvprd(n, x, incx);
}
inline void blasDgeva(char* trans, BlasInt* n, BlasInt* m, double* alpha,
					double* a, double* b) {
	dgeva(trans, n, m, alpha, a, b);
}
inline void blasDgevm(char* trans, BlasInt* m, BlasInt* n, double* alpha,
					double* a, double* b) {
	dgevm(trans, m, n, alpha, a, b);
}
inline void blasDgesu(char* trans, BlasInt* n, BlasInt* m, double* alpha,
					double* a, double* beta, double* b) {
	dgesu(trans, n, m, alpha, a, beta, b);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* BLAS_EXT_H_ */
