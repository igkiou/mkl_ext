/*
 * mkl_ext//mkl_ext/mexfiles/dchrk.cpp/dchrk.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: igkiou
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "blas_ext.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Check number of input arguments */
	if (nrhs != 5) {
		mexErrMsgTxt("Five input arguments are required.");
    }

	/* Check number of output arguments */
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
    }

	double *C = (double *) mxGetData(prhs[0]);
	BlasInt isU = (BlasInt) * (double *) mxGetData(prhs[1]);
	double *A = (double *) mxGetData(prhs[2]);
	BlasInt transN = (BlasInt) * (double *) mxGetData(prhs[3]);
	double *a = (double *) mxGetData(prhs[4]);

	BlasInt N;
	BlasInt K;
	if (transN == 1) {
		N = (BlasInt) mxGetM(prhs[2]);
		K = (BlasInt) mxGetN(prhs[2]);
	} else {
		N = (BlasInt) mxGetN(prhs[2]);
		K = (BlasInt) mxGetM(prhs[2]);
	}

	if (mxGetNumberOfElements(prhs[4]) != 1) {
		mexErrMsgTxt("a must be a scalar.");
	}
	if (mxGetM(prhs[0]) != N) {
		mexErrMsgTxt("A must be a N x N matrix.");
	}
	if (mxGetM(prhs[0]) != mxGetN(prhs[0])) {
		mexErrMsgTxt("A must be a square matrix.");
	}

	plhs[0] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
	double *L = (double *) mxGetData(plhs[0]);
	char uplo;
	if (isU == 1) {
		 uplo = 'U';
	} else {
		uplo = 'L';
	}
	char trans;
	if (transN == 1) {
		trans = 'N';
	} else {
		trans = 'T';
	}

	BlasInt ione = 1;
	double *work = (double *) malloc(N * sizeof(double));
	BlasInt info = 0;
	BlasInt Nsq = N * N;
	dcopy(&Nsq, C, &ione, L, &ione);
	dchrk(&uplo, &trans, &N, &K, a, A, L, work, &info);

	plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double *temp = (double *) mxGetData(plhs[1]);
	*temp = (double) info;
	free(work);
}
