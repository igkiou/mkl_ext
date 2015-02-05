/*
 * linalg_benchmarking//linalg_benchmarking/armadillo-3.900.6/examples/gemm_arma_test.cpp/gemm_arma_test.cpp
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

	arma::mat Amat(A, M, N, false, true);
	arma::mat Bmat(B, O, P, false, true);
	arma::mat Cmat(C, Mc, Nc, false, true);

	for (int iter = 0; iter < numIters; ++iter) {
		if ((transA == 0) && (transB == 0)) {
			Cmat = Amat * Bmat;
		} else if ((transA == 0) && (transB == 1)) {
			Cmat = Amat * Bmat.t();
		} else if ((transA == 1) && (transB == 0)) {
			Cmat = Amat.t() * Bmat;
		} else if ((transA == 1) && (transB == 1)) {
			Cmat = Amat.t() * Bmat.t();
		}
	}
}
