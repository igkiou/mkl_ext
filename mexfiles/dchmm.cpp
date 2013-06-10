/*
 * mkl_ext//mkl_ext/mexfiles/dchmm.cpp/dchmm.cpp
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

	double *A = (double *) mxGetData(prhs[0]);
	BlasInt isU = (BlasInt) * (double *) mxGetData(prhs[1]);
	double *B = (double *) mxGetData(prhs[2]);
	BlasInt sideL = (BlasInt) * (double *) mxGetData(prhs[3]);
	double *a = (double *) mxGetData(prhs[4]);

	BlasInt M = (BlasInt) mxGetM(prhs[2]);
	BlasInt N = (BlasInt) mxGetN(prhs[2]);
	BlasInt MN = M * N;

	if (mxGetNumberOfElements(prhs[4]) != 1) {
		mexErrMsgTxt("a must be a scalar.");
	}
	if (sideL == 1) {
		if (mxGetM(prhs[0]) != M) {
			mexErrMsgTxt("A must be a M x M matrix.");
		}
	} else {
		if (mxGetM(prhs[0]) != N) {
			mexErrMsgTxt("A must be a N x N matrix.");
		}
	}
	if (mxGetM(prhs[0]) != mxGetN(prhs[0])) {
		mexErrMsgTxt("A must be a square matrix.");
	}

	plhs[0] = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
	double *R = (double *) mxGetData(plhs[0]);
	BlasInt ione = 1;
	dcopy(&MN, B, &ione, R, &ione);
	char uplo;
	if (isU == 1) {
		uplo = 'U';
	} else {
		uplo = 'L';
	}
	char side;
	if (sideL == 1) {
		side = 'L';
	} else {
		side = 'R';
	}
	dchmm(&side, &uplo, &M, &N, a, A, R);
}



