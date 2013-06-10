/*
 * mkl_ext//mkl_ext/mexfiles/dchmv.cpp/dchmv.cpp
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
	if (nrhs != 3) {
		mexErrMsgTxt("Three input arguments are required.");
    }

	/* Check number of output arguments */
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
    }

	double *A = (double *) mxGetData(prhs[0]);
	BlasInt isU = (BlasInt) * (double *) mxGetData(prhs[1]);
	double *v = (double *) mxGetData(prhs[2]);

	BlasInt N = (BlasInt) mxGetM(prhs[0]);

	if (mxGetN(prhs[0]) != N) {
		mexErrMsgTxt("A must be a square matrix.");
	} else if (mxGetM(prhs[2]) != N) {
		mexErrMsgTxt("v must be of the same length as a side of A.");
	} else if (mxGetN(prhs[2]) != 1) {
		mexErrMsgTxt("v must be a vector.");
	}

	plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
	double *y = (double *) mxGetData(plhs[0]);
	char uplo;
	if (isU == 1) {
		 uplo = 'U';
	} else {
		uplo = 'L';
	}
	BlasInt ione = 1;
	dcopy(&N, v, &ione, y, &ione);
	dchmv(&uplo, &N, A, y, &ione);
}
