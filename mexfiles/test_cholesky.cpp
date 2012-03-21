/*
 * l1_distance_mex.c
 *
 *  Created on: Mar 16, 2012
 *      Author: igkiou
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mkl_cholesky.h"
#include "mex.h"
#include "matrix.h"

void print_matrix(double *matrix, MKL_INT M, MKL_INT N) {

	MKL_INT iterM, iterN;
	for (iterM = 0; iterM < M; ++iterM) {
		for (iterN = 0; iterN < N; ++iterN) {
			mexPrintf("%g ", matrix[iterN * M + iterM]);
		}
		mexPrintf("\n");
	}
}

void print_matrix_int(MKL_INT *matrix, MKL_INT M, MKL_INT N) {

	MKL_INT iterM, iterN;
	for (iterM = 0; iterM < M; ++iterM) {
		for (iterN = 0; iterN < N; ++iterN) {
			mexPrintf("%d ", matrix[iterN * M + iterM]);
		}
		mexPrintf("\n");
	}
}

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
	double *v = (double *) mxGetData(prhs[1]);
	double *a = (double *) mxGetData(prhs[2]);

	MKL_INT N = (MKL_INT) mxGetM(prhs[0]);

	if (mxGetN(prhs[0]) != N) {
		mexErrMsgTxt("A must be a square matrix.");
	} else if (mxGetM(prhs[1]) != N) {
		mexErrMsgTxt("v must be of the same length as a side of A.");
	} else if (mxGetN(prhs[1]) != 1) {
		mexErrMsgTxt("v must be a vector.");
	}

	plhs[0] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
	double *L = (double *) mxGetData(plhs[0]);
	char uplo = 'U';
	mexPrintf("gkiou1\n");
	MKL_INT incx = 1;
	MKL_INT lda = N;
	double *work = (double *) malloc(3 * N * sizeof(double));
	MKL_INT lwork = 3 * N;
	MKL_INT info;
	MKL_INT Nsq = N * N;
	dcopy(&Nsq, A, &incx, L, &incx);
	dchr(&uplo, &N, a, v, &incx, L, &lda, work, &lwork, &info);

	plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double *temp = (double *) mxGetData(plhs[1]);
	*temp = (int) info;
	free(work);
}
