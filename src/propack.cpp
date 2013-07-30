/*
 * mkl_ext//mkl_ext/src/propack.cpp/propack.cpp
 *
 *  Created on: Jun 17, 2013
 *      Author: igkiou
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "propack.h"

void dlansvd(char* jobu, char* jobvt, BlasInt *m, BlasInt *n, double *a,
			double *u, double *vt, BlasInt *k, BlasInt *job);

void dlanbpro(double* a, BlasInt* m, BlasInt* n, BlasInt* k, double* r0,
			double* u, double* b, double* vt, BlasInt* f, double* work,
			BlasInt* ierr) {

	const BlasInt M = *m;
	const BlasInt N = *n;
	const BlasInt K = *k;
	const BlasInt F = *f;

	const double m2 = 1.5;
	const double n2 = 1.5;
	const char lamchArg = 'E';
	const double eps = dlamch(&lamchArg);
	const double sqrtEps = sqrt(eps);
	const double delta = sqrtEps / sqrt((double) k);
	const double eta = delta * sqrt(sqrtEps);
	const bool cgs = false;
	const bool elr = true;
	const double gamma = 1 / sqrt(2.0);
	const BlasInt onesided = 0;
	double t = 0;

	double FUDGE = 1.01;
	BlasInt npu = 0;
	BlasInt npv = 0;
	*ierr = 0;

	double* p = r0;
	double* alpha = &work[0];
	double* beta = &work[K + 1];
	/*
	 * TODO: Check definition of r.
	 */
	double* r = &work[2 * K + 2];
	const BlasInt ione = 1;
	double anorm;
	bool est_anorm;
	bool force_reorth;

	double* mu;
	double* nu;
	double* numax;

	/*
	 * TODO: Add preparation for Lanczos iteration here.
	 */

	for (BlasInt lanIter = F + 1; lanIter <= K; ++lanIter) {

		dcopy(m, p, &ione, &u[(lanIter - 1) * m], &ione);
		if (beta[lanIter] != 0) {
			const double scale = 1 / beta[lanIter - 1];
			dscal(m, &scale, &u[(lanIter - 1) * m], &ione);
		}

		if (lanIter == 6) {
			const BlasInt length = lanIter - 1;
			const BlasInt incx = 	lanIter + 1;
			dcopy(&length, &alpha[0], &ione, &b[0], &incx);
			dcopy(&length, &beta[1], &ione, &b[1], &incx);
			/*
			 * TODO: Fill in 2-norm computation here.
			 */
//			anorm = FUDGE *
			est_anorm = false;
		}

		if (lanIter == 1) {
			const char trans = 'T';
			const double alphaVal = 1.0;
			const double betaVal = 0;
			dgemv(&trans, &M, &N, &alphaVal, a, &M, u, &ione, &betaVal, r,
				&ione);
			alpha[0] = dnrm2(&N, r, &ione);
			if (est_anorm) {
				anorm = FUDGE * alpha[0];
			}
		} else {
			dcopy(&N, &vt[lanIter - 2], &K, r, &ione);
			const char trans = 'T';
			const double alphaVal = 1.0;
			const double betaVal = -1.0;
			dgemv(&trans, &M, &N, &alphaVal, a, &M, u, &ione, &betaVal, r,
				&ione);
			alpha[lanIter - 1] = dnrm2(&N, r, &ione);

			if ((alpha[lanIter - 1] < gamma * beta[lanIter - 1]) && (elr)) {
				double normold = alpha[lanIter - 1];
				bool stop = false;
				while (!stop) {
					const BlasInt incy = M;
					t = ddot(&N, r, &ione, &vt[lanIter - 2], &incy);
					const double scale = -t;
					daxpy(&N, &scale, &vt[lanIter - 2], &incy, r, &ione);
					alpha[lanIter - 1] = dnrm2(&N, r, &ione);
					if (beta[lanIter - 1] != 0) {
						beta[lanIter - 1] += t;
					}
					if (alpha[lanIter - 1] >= gamma * normold) {
						stop = true;
					} else {
						normold = alpha[lanIter - 1];
					}
				}
			}
			if (est_anorm) {
				if (lanIter == 2) {
					double anormTemp = FUDGE * sqrt(alpha[0] * alpha[0]
					                               + beta[1] * beta[1]
					                               + alpha[1] * beta[1]);
					if (anormTemp > anorm) {
						anorm = anormTemp;
					}
				} else {
					double anormTemp = FUDGE
								* sqrt(alpha[lanIter - 2] * alpha[lanIter - 2]
								+ beta[lanIter - 1] * beta[lanIter - 1]
								+ alpha[lanIter - 2] * beta[lanIter - 2]
								+ alpha[lanIter - 1] * beta[lanIter - 1]);
					if (anormTemp > anorm) {
						anorm = anormTemp;
					}
				}
			}

			if (alpha[lanIter - 1] != 0) {
				update_nu(nu, mu, lanIter, alpha, beta, anorm);
				BlasInt length = lanIter - 1;
				BlasInt maxIndex = idamax(&length, nu, &ione);
				numax[lanIter] = nu[maxIndex - 1];
			}

			if (elr) {
				nu[lanIter - 2] = n2 * eps;
			}

			if ((onesided != -1)
					&& ((numax[lanIter - 1] > delta) || force_reorth)
					&& alpha[lanIter - 1] != 0) {
				if (eta == 0) {
					
				}
			}

		}
	}
}

void update_nu(double* nu, double* mu, const BlasInt lanIter,
			const double* alpha, const double* beta, const double anorm) {
	const double ainv = 1 / alpha[lanIter - 1];
	const char lamchArg = 'E';
	const double eps = dlamch(&lamchArg);
	const double eps1 = 100.0 * eps / 2.0;
	if (lanIter > 1) {
		for (BlasInt iterK = 1; iterK <= lanIter - 1; ++iterK) {
			double T = eps1 * (hypot(alpha[iterK - 1], beta[iterK])
							+ hypot(alpha[lanIter - 1], beta[lanIter - 1]));
			T += eps1 * anorm;
			nu[iterK - 1] = beta[iterK] * mu[iterK]
			              + alpha[iterK - 1] * mu[iterK - 1]
			              + beta[lanIter - 1] * nu[iterK - 1];
			nu[iterK - 1] = ainv
						* (nu[iterK - 1] + copysign(1.0, nu[iterK - 1]) * T);
		}
	}
	nu[lanIter - 1] = 1.0;
}


//[U,B,V,p,ierr,w] = lanbpro_modified_nocomplex(A,j,p,options,U,B,V,anorm);
//[U_k,B_k,V_k,R,ierr,work] = LANBPRO(A,K,R0,OPTIONS,U_old,B_old,V_old)



void dmgs(BlasInt* n, BlasInt* k, double* V, BlasInt* ldv, double* vnew,
		BlasInt* index) {

	const BlasInt LDV = *ldv;
	for (BlasInt iterI = 0; iterI < *k; ++iterI) {
		BlasInt idx = index[iterI] - 1;  /* -1 because MATLAB uses 1-based indices */
		double s = 0.0;

		for (BlasInt iterJ = 0; iterJ < *n; ++iterJ) {
			/*s += V[j,idx]*vnew[j]; */ 
			s += V[idx *LDV + iterJ] * vnew[iterJ];
		}

		for (BlasInt iterJ = 0; iterJ < *n; ++iterJ) {
			/* vnew[j] -= s*V[iterJ,idx]; */ /* Fortran is row-indexed */
			vnew[iterJ] -= s * V[idx * LDV + iterJ];
		}
	}
}

void reorth(BlasInt* n, BlasInt* k, double* V, BlasInt* ldv, double* vnew, 
		double* normv, BlasInt* index, double* alpha, double* work,
		BlasInt* iflag, BlasInt* nre) {

	const BlasInt ione = 1;
	const BlasInt N = *n;
	const BlasInt K = *k;
	const BlasInt LDV = *ldv;
	const char transpose = 'T';
	const char normal = 'N';
	double normv_old;
	const BlasInt maxTry = 4;

	const double done = 1.0;
	const double dmone = -1.0;
	const double dzero = 0.0;

//	BlasInt workflag = 0;
//	if (work == NULL) {
//		work = mxMalloc(* k * sizeof(double));
//		workflag = 1;
//	}
	
	/* Hack: if index != 1:k, we do MGS to avoid reshuffling */
	if (*iflag == 1) {
		for (BlasInt iter = 0; iter < *k; ++iter){
			if (index[iter] != (iter + 1)) {
				*iflag = 0;
				break;
			}
		}
	}
	normv_old = 0;
	*nre = 0;  
	*normv = 0.0;

	while ((*normv < *alpha * normv_old) || (*nre == 0)) {
		if (*iflag == 1) {
			dgemv(&transpose, &N, &K, &done, V, &LDV, vnew, &ione, &dzero, work, &ione);
			dgemv(&normal, &N, &K, &dmone, V, &LDV, work, &ione, &done, vnew, &ione);
		} else {
			/* MGS */
			dmgs(n, k, V, ldv, vnew, index);
		}
		normv_old = *normv; 
		/* following line works! */
		*normv = dnrm2(&N, vnew, &ione);

		*nre = *nre + 1;

		if (*nre > maxTry) {
			/* vnew is numerically in span(V) --> return vnew as all zeros */
			*normv = 0.0;
			for (BlasInt iter = 0; iter< *n ; ++iter)
				vnew[iter] = 0.0;
			return;
		}
	}
}


