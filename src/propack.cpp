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
	double* anorm;
	BlasInt est_anorm;

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
			est_anorm = 0;
		}

		if (lanIter == 1) {
			const char trans = 'T';
			const double alphaVal = 1.0;
			const double betaVal = 0;
			dgemv(&trans, &M, &N, &alphaVal, a, &M, u, &ione, &betaVal, r, &ione);
			alpha[0] = dnrm2(&N, r, &ione);
			if (est_anorm > 0) {
				anorm = FUDGE * alpha[0];
			}
		} else {
			dcopy(&N, &vt[lanIter - 2], &K, r, &ione);
			const char trans = 'T';
			const double alphaVal = 1.0;
			const double betaVal = -1.0;
			dgemv(&trans, &M, &N, &alphaVal, a, &M, u, &ione, &betaVal, r, &ione);
			alpha[lanIter - 1] = dnrm2(&N, r, &ione);

			if ((alpha[lanIter - 1] < gamma * beta[lanIter - 1]) && (elr)) {
				double normOld = alpha[lanIter - 1];
			}
		}
	}
}

//[U,B,V,p,ierr,w] = lanbpro_modified_nocomplex(A,j,p,options,U,B,V,anorm);
//[U_k,B_k,V_k,R,ierr,work] = LANBPRO(A,K,R0,OPTIONS,U_old,B_old,V_old)



