/*
 * linalg_benchmarking//linalg_benchmarking/benchmarking_examples/symm_mkl_test.cpp/symm_mkl_test.cpp
 *
 *  Created on: Jul 10, 2013
 *      Author: igkiou
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "armadillo"
#include "mex.h"
#include "matrix.h"
#include "mkl.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Check number of input arguments */
	if (nrhs != 5) {
		mexErrMsgTxt("Five input arguments are required.");
    }

	/* Check number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
    }

	double *A = (double *) mxGetData(prhs[0]);
	double *B = (double *) mxGetData(prhs[1]);
	const int sideL = (int) * (double *) mxGetData(prhs[2]);
	const int uploU = (int) * (double *) mxGetData(prhs[3]);
	const int numIters = (int) * (double *) mxGetData(prhs[4]);

	const int M = (int) mxGetM(prhs[0]);
	const int N = (int) mxGetM(prhs[1]);
	const int O = (int) mxGetN(prhs[1]);
	int Mc;
	int Nc;

	if (M != (int) mxGetM(prhs[0])) {
		mexErrMsgTxt("First (symmetric) matrix must be square.");
	}

	if (sideL == 1) {
		if (N != M) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = M;
		Nc = O;
	} else if (sideL == 0) {
		if (O != M) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = N;
		Nc = M;
	}
	plhs[0] = mxCreateNumericMatrix(Mc, Nc, mxDOUBLE_CLASS, mxREAL);
	double *C = (double *) mxGetData(plhs[0]);
	const char sideLc = (sideL == 1)?('L'):('R');
	const char uploUc = (uploU == 1)?('U'):('L');
	const double alpha = 1.0;
	const double beta = 0;

	for (int iter = 0; iter < numIters; ++iter) {
		dsymm(&sideLc, &uploUc, &Mc, &Nc, &alpha, A, &M, B, &Mc, &beta, C, &Mc);
	}
}



