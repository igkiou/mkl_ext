/*
 * mkl_ext//mkl_ext/linalg_benchmarking/benchmarking_examples/gemm_blaze_test.cpp/gemm_blaze_test.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: igkiou
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "blaze/Blaze.h"
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
	const int transA = (int) * (double *) mxGetData(prhs[2]);
	const int transB = (int) * (double *) mxGetData(prhs[3]);
	const int numIters = (int) * (double *) mxGetData(prhs[4]);

	const int M = (int) mxGetM(prhs[0]);
	const int N = (int) mxGetN(prhs[0]);
	const int O = (int) mxGetM(prhs[1]);
	const int P = (int) mxGetN(prhs[1]);
	int Mc;
	int Nc;
	int Kc;

	if ((transA == 0) && (transB == 0)) {
		if (N != O) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = M;
		Nc = P;
		Kc = N;
	} else if ((transA == 0) && (transB == 1)) {
		if (N != P) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = M;
		Nc = O;
		Kc = N;
	} else if ((transA == 1) && (transB == 0)) {
		if (M != O) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = N;
		Nc = P;
		Kc = M;
	} else if ((transA == 1) && (transB == 1)) {
		if (M != P) {
			mexErrMsgTxt("Incompatible matrix dimensions.");
		}
		Mc = N;
		Nc = O;
		Kc = M;
	}
	plhs[0] = mxCreateNumericMatrix(Mc, Nc, mxDOUBLE_CLASS, mxREAL);
	double *C = (double *) mxGetData(plhs[0]);

	blaze::DynamicMatrix<double, blaze::columnMajor> Amat(M, N, A);
	blaze::DynamicMatrix<double, blaze::columnMajor> Bmat(O, P, B);
	blaze::DynamicMatrix<double, blaze::columnMajor> Cmat(Mc, Nc, C);

	for (int iter = 0; iter < numIters; ++iter) {
		if ((transA == 0) && (transB == 0)) {
			Cmat = Amat * Bmat;
		} else if ((transA == 0) && (transB == 1)) {
			Cmat = Amat * blaze::trans(Bmat);
		} else if ((transA == 1) && (transB == 0)) {
			Cmat = blaze::trans(Amat) * Bmat;
		} else if ((transA == 1) && (transB == 1)) {
			Cmat = blaze::trans(Amat) *blaze::trans(Bmat);
		}
	}
}
