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

typedef struct SFMTEngineState {
	union {
		/** the 128-bit internal state array */
		w128_t sfmt[N];
		/** the 32bit integer version of the 128-bit internal state array */
		uint32_t psfmt32[N32];
		/** the 64bit integer version of the 128-bit internal state array */
		uint64_t psfmt64[N64];
	};

	/** index counter to the 32-bit internal state array */
	int idx;

	/** a parity check vector which certificate the period of 2^{MEXP} */
	const static uint32_t parity[4];

	/** Hash of the SFMT parameters */
	const static uint32_t s_magic;

	/* Default constructor, set the index to an invalid value */
	State() : idx(-1) {}
} SFMTEngineState;

inline bool isInitialized(const SFMTEngineState* state) {
	return state->idx >= 0;
}

void init_gen_rand(SFMTEngineState* state, uint64_t seed);

void init_by_array(SFMTEngineState* state, const uint32_t *init_key,
					int key_length);

inline uint64_t gen_rand64(SFMTEngineState* state) {
	if (state->idx >= N32) {
		gen_rand_all(state);
		state->idx = 0;
	}

	uint64_t r = state->psfmt64[state->idx / 2];
	state->idx += 2;
	return r;
}

static inline uint32_t func1(uint32_t x) {
	return (x ^ (x >> 27)) * (uint32_t)1664525UL;
}

static inline uint32_t func2(uint32_t x) {
	return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
}

void period_certification(SFMTEngineState* state) {
	int inner = 0;
	int i, j;
	uint32_t work;

	for (i = 0; i < 4; ++i)
		inner ^= state->psfmt32[i] & state->parity[i];
	for (i = 16; i > 0; i >>= 1)
		inner ^= inner >> i;
	inner &= 1;
	/* check OK */
	if (inner == 1) {
		return;
	}
	/* check NG, and modification */
	for (i = 0; i < 4; ++i) {
		work = 1;
		for (j = 0; j < 32; ++j) {
			if ((work & state->parity[i]) != 0) {
				state->psfmt32[i] ^= work;
				return;
			}
			work = work << 1;
		}
	}
}

inline void gen_rand_all(SFMTEngineState* state) {
	int i;
	__m128i r, r1, r2, mask;
	mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

	r1 = _mm_load_si128(&(state->sfmt[N - 2].si));
	r2 = _mm_load_si128(&(state->sfmt[N - 1].si));
	for (i = 0; i < N - POS1; ++i) {
		r = mm_recursion(state->sfmt[i].si, state->sfmt[i + POS1].si, r1, r2, mask);
		_mm_store_si128(&(state->sfmt[i].si), r);
		r1 = r2;
		r2 = r;
	}
	for (; i < N; ++i) {
		r = mm_recursion(&(state->sfmt[i].si), state->sfmt[i + POS1 - N].si, r1, r2, mask);
		_mm_store_si128(&(state->sfmt[i].si), r);
		r1 = r2;
		r2 = r;
	}
}



















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
