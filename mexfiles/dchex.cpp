/*
 * mkl_ext//mkl_ext/mexfiles/dchex.cpp/dchex.cpp
 *
 *  Created on: Jun 16, 2013
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
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
    }

	double *A = (double *) mxGetData(prhs[0]);
	BlasInt k = (BlasInt) * (double *) mxGetData(prhs[1]);
	BlasInt l = (BlasInt) * (double *) mxGetData(prhs[2]);
	BlasInt job = (BlasInt) * (double *) mxGetData(prhs[3]);

	BlasInt N = (BlasInt) mxGetM(prhs[0]);

	if (mxGetN(prhs[0]) != N) {
		mexErrMsgTxt("A must be a square matrix.");
	} else if ((k < 1) || (k > N)) {
		mexErrMsgTxt("Incorrect value for K.");
	} else if ((l < k) || (l > N)) {
		mexErrMsgTxt("Incorrect value for L.");
	}

	plhs[0] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
	double *L = (double *) mxGetData(plhs[0]);

	double *work = (double *) malloc(2 * N * sizeof(double));
	BlasInt ione = 1;
	BlasInt Nsq = N * N;
	dcopy(&Nsq, A, &ione, L, &ione);
	dchex(L, &N, &k, &l, work, &job);
	free(work);
}
