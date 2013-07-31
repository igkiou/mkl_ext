/*
 * mkl_ext//mkl_ext/mexfiles/dchr.cpp/dchr.cpp
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
	if (nrhs != 4) {
		mexErrMsgTxt("Four input arguments are required.");
    }

	/* Check number of output arguments */
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
    }

	double *A = (double *) mxGetData(prhs[0]);
	BlasInt isU = (BlasInt) * (double *) mxGetData(prhs[1]);
	double *v = (double *) mxGetData(prhs[2]);
	double *a = (double *) mxGetData(prhs[3]);

	BlasInt N = (BlasInt) mxGetM(prhs[0]);

	if (mxGetN(prhs[0]) != N) {
		mexErrMsgTxt("A must be a square matrix.");
	} else if (mxGetM(prhs[2]) != N) {
		mexErrMsgTxt("v must be of the same length as a side of A.");
	} else if (mxGetN(prhs[2]) != 1) {
		mexErrMsgTxt("v must be a vector.");
	} else if (mxGetNumberOfElements(prhs[3]) != 1) {
		mexErrMsgTxt("a must be a scalar.");
	}

	plhs[0] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
	double *L = (double *) mxGetData(plhs[0]);
	char uplo;
	if (isU == 1) {
		 uplo = 'U';
	} else {
		uplo = 'L';
	}
	BlasInt ione = 1;
	double *work = (double *) malloc(N * sizeof(double));
	BlasInt info = 0;
	BlasInt Nsq = N * N;
	dcopy(&Nsq, A, &ione, L, &ione);
	dchr(&uplo, &N, a, v, &ione, L, work, &info);

	plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double *temp = (double *) mxGetData(plhs[1]);
	*temp = (double) info;
	free(work);
}
