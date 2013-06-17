/*
 * mkl_ext//mkl_ext/include/sfmt_header.h/sfmt_header.h
 *
 *  Created on: Jun 17, 2013
 *      Author: igkiou
 */

#ifndef SFMT_HEADER_H_
#define SFMT_HEADER_H_

#include <stdint.h>
#include "blas_header.h"

#define Epsilon 1e-7
#define ShadowEpsilon 1e-5
#define DeltaEpsilon 1e-3f

/* Assumed L1 cache line size for alignment purposes */
#if !defined(L1_CACHE_LINE_SIZE)
#define L1_CACHE_LINE_SIZE 64
#endif

#ifdef M_E
#undef M_E
#endif

#ifdef M_PI
#undef M_PI
#endif

#ifdef INFINITY
#undef INFINITY
#endif

#define ONE_MINUS_EPS_FLT 0x1.fffffep-1f
#define ONE_MINUS_EPS_DBL 0x1.fffffffffffff7p-1

#define M_E           2.71828182845904523536
#define M_PI          3.14159265358979323846
#define INV_PI        0.31830988618379067154
#define INV_TWOPI     0.15915494309189533577
#define INV_FOURPI    0.07957747154594766788
#define SQRT_TWO      1.41421356237309504880
#define INV_SQRT_TWO  0.70710678118654752440
#define ONE_MINUS_EPS ONE_MINUS_EPS_DBL

typedef uint64_t SFMTEngineSeedType;

class SSEEngine {
public:
	SSEEngine();

	explicit SSEEngine(const IndexType seedValue);

	void seed(const BlasInt seedValue);




	inline double operator()() {
		/* Trick from MTGP: generate an uniformly distributed
		   single precision number in [1,2) and subtract 1. */
		union {
			uint64_t u;
			double d;
		} x;
		x.u = (nextULong() >> 12) | 0x3ff0000000000000ULL;
		return x.d - 1.0;
	}

	~SSEEngine();

private:
	/// Return an integer on the [0, 2^63-1]-interval
	uint64_t nextULong();

	/// Return an integer on the [0, n)-interval
	uint32_t nextUInt(uint32_t n);

	/// Return an integer on the [0, n)-interval
	size_t nextSize(size_t n);

	/// Return a floating point value on the [0, 1) interval
	double nextdouble();

	/// Return a normally distributed value
	double nextStandardNormal();

	struct State;
	State *mt;
};

#endif /* SFMT_HEADER_H_ */
