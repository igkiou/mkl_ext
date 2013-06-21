/*
 * mkl_ext//mkl_ext/src/rng.cpp/rng.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: igkiou
 */

#include <math.h>
#include <stdint.h>

#include "rng.h"
#define UNUSED(x) (void)(x)
#define M_PI 3.14159265358979323846

#if defined(USE_RNG_VSL)

void rngInit(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info) {
	*info = vslNewStream(&(rng->m_stream), kVSLBRNGMethod, *seedValue);
}

void drunif(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {
	UNUSED(isAligned);
	*info = (*n > 0)?(vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rng->m_stream,
								*n, buffer, 0.0, 1.0))
					:(0);
}

void drnorm(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {
	UNUSED(isAligned);
#if defined(USE_RNG_BOX_MULLER)
	if (*n == 1) {
		*info = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, rng->m_stream,
							1, buffer, 0.0, 1.0);
	} else if (*n > 1) {
		*info = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, rng->m_stream,
							*n, buffer, 0.0, 1.0);
	} else {
		*info = 0;
	}
#elif defined(USE_RNG_MARSAGLIA)
	if (*n > 0) {
		for (BlasInt iter = 0; iter < *n; ++iter) {
			double x[2];
			double r;
			do {
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rng->m_stream, 2, x,
							0.0, 1.0);
				x[0] = 2.0 * x[0] - 1.0;
				x[1] = 2.0 * x[1] - 1.0;
				r = x[0] * x[0] + x[1] * x[1];
			} while (r >= 1 || r == 0);
			buffer[iter] = x[0] * sqrt(- 2.0 * log(r) / r);
		}
	} else {
		*info = 0;
	}
#endif
}

void irunif(RngEngine* rng, BlasInt* r, BlasInt* buffer, BlasInt* n,
			BlasInt* isAligned, RngErrorType* info) {
	UNUSED(isAligned);
	*info = (*n > 0)?(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rng->m_stream,
								*n, buffer, 0, *r))
					:(0);
}

void rngDestroy(RngEngine* rng, RngErrorType* info) {
	*info = (rng->m_stream != NULL)?(vslDeleteStream(&(rng->m_stream))):(0);
}
#elif defined(USE_RNG_SFMT)

void rngInit(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info) {
	sfmt_init_gen_rand(&(rng->m_sfmt), *seedValue);
	*info = 0;
}

void drunif(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {
	UNUSED(isAligned);
	if (*n > 0) {
		for (BlasInt iter = 0; iter < *n; ++iter) {
			buffer[iter] = sfmt_genrand_real2(&(rng->m_sfmt));
		}
	}
	*info = 0;
}

void drnorm(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {
	UNUSED(isAligned);
	if (*n > 0) {
#if defined(USE_RNG_BOX_MULLER)
		for (BlasInt iter = 0; iter < *n; iter += 2) {
			double x = 2.0 * M_PI * sfmt_genrand_real2(&(rng->m_sfmt));
			double z = sqrt(-log(sfmt_genrand_real2(&(rng->m_sfmt))));
			double c = cos(x);
			buffer[iter] = z * c;
			if (iter + 1 < *n) {
				double s = sin(x);
				buffer[iter + 1] = z * s;
			}
		}
#elif defined(USE_RNG_MARSAGLIA)
		for (BlasInt iter = 0; iter < *n; ++iter) {
			double x, y, r;
			do {
				x = 2.0 * sfmt_genrand_real2(&(rng->m_sfmt)) - 1.0;
				y = 2.0 * sfmt_genrand_real2(&(rng->m_sfmt)) - 1.0;
				r = x * x + y * y;
			} while (r >= 1 || r == 0);
			buffer[iter] = x * sqrt(- 2.0 * log(r) / r);
		}
#endif
	}
	*info = 0;
}

void irunif(RngEngine* rng, BlasInt* r, BlasInt* buffer, BlasInt* n,
			BlasInt* isAligned, RngErrorType* info) {
	UNUSED(isAligned);
	if (*n > 0) {
		for (BlasInt iter = 0; iter < *n; ++iter) {
			buffer[iter] = (BlasInt) sfmt_genrand_uint32(&(rng->m_sfmt));
		}
	}
	*info = 0;
}

void rngDestroy(RngEngine* rng, RngErrorType* info) {
	UNUSED(rng);
	*info = 0;
}
#elif defined(USE_RNG_DSFMT)

void rngInit(RngEngine* rng, RngSeedType* seedValue, RngErrorType* info) {
	dsfmt_init_gen_rand(&(rng->m_dsfmt), *seedValue);
	*info = 0;
}

void drunif(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {

	if ((*isAligned == 0) || (*n == 1)) {
		for (BlasInt iter = 0; iter < *n; ++iter) {
			buffer[iter] = dsfmt_genrand_close_open(&(rng->m_dsfmt));
		}
	} else if (*n > 1) {
		dsfmt_fill_array_close_open(&(rng->m_dsfmt), buffer, *n);
	}
	*info = 0;
}

void drnorm(RngEngine* rng, double* buffer, BlasInt* n, BlasInt* isAligned,
			RngErrorType* info) {
	UNUSED(isAligned);
	if (*n > 0) {
#if defined(USE_RNG_BOX_MULLER)
		for (BlasInt iter = 0; iter < *n; iter += 2) {
			double x = 2.0 * M_PI * dsfmt_genrand_close_open(&(rng->m_dsfmt));
			double z = sqrt(-log(dsfmt_genrand_close_open(&(rng->m_dsfmt))));
			double c = cos(x);
			buffer[iter] = z * c;
			if (iter + 1 < *n) {
				double s = sin(x);
				buffer[iter + 1] = z * s;
			}
		}
#elif defined(USE_RNG_MARSAGLIA)
		for (BlasInt iter = 0; iter < *n; ++iter) {
			double x, y, r;
			do {
				x = 2.0 * dsfmt_genrand_close_open(&(rng->m_dsfmt)) - 1.0;
				y = 2.0 * dsfmt_genrand_close_open(&(rng->m_dsfmt)) - 1.0;
				r = x * x + y * y;
			} while (r >= 1 || r == 0);
			buffer[iter] = x * sqrt(- 2.0 * log(r) / r);
		}
#endif
	}
	*info = 0;
}

void irunif(RngEngine* rng, BlasInt* r, BlasInt* buffer, BlasInt* n,
			BlasInt* isAligned, RngErrorType* info) {
	UNUSED(isAligned);
	if (*n > 0) {
		for (BlasInt iter = 0; iter < *n; ++iter) {
			double x = dsfmt_genrand_close_open(&(rng->m_dsfmt));
			buffer[iter] = (BlasInt) floor(*r * x);
		}
	}
	*info = 0;
}

void rngDestroy(RngEngine* rng, RngErrorType* info) {
	UNUSED(rng);
	*info = 0;
}
#endif
