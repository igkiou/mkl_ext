/*
 * mkl_ext//mkl_ext/mexfiles/drunif.cpp/drunif.cpp
 *
 *  Created on: Jun 21, 2013
 *      Author: igkiou
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "blas_ext.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Check number of input arguments */
	if (nrhs != 2) {
		mexErrMsgTxt("Two input arguments are required.");
	}

	/* Check number of output arguments */
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
	}

	BlasInt M = (BlasInt) * (double *) mxGetData(prhs[0]);
	BlasInt N = (BlasInt) * (double *) mxGetData(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
	double *data = (double *) mxGetData(plhs[0]);

	RngEngine rng;
	RngErrorType info;
	RngSeedType seedValue = (RngSeedType) time(0);
	rngInit(&rng, &seedValue, &info);
	BlasInt MN = M * N;
	BlasInt isAligned = 0;
	drunif(&rng, data, &MN, &isAligned, &info);
	rngDestroy(&rng, &info);
}
