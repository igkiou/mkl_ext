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
			double* u, double* b, double* vt, double* work, BlasInt* ierr) {

	BlasInt M = *m;
	BlasInt N = *n;
	BlasInt K = *k;

	double m2 = 1.5;
	double n2 = 1.5;
	char lamchArg = 'E';
	double eps = dlamch(&lamchArg);
	double sqrtEps = sqrt(eps);
	double delta = sqrtEps / sqrt((double) k);
	double eta = delta * sqrt(sqrtEps);
	BlasInt cgs = 0;
	BlasInt elr = 2;
	double gamma = 1 / sqrt(2.0);
	BlasInt onesided = 0;

	double FUDGE = 1.01;
	BlasInt npu = 0;
	BlasInt npv = 0;
	*ierr = 0;

}

//[U,B,V,p,ierr,w] = lanbpro_modified_nocomplex(A,j,p,options,U,B,V,anorm);
//[U_k,B_k,V_k,R,ierr,work] = LANBPRO(A,K,R0,OPTIONS,U_old,B_old,V_old)



