/* Small helpers of the embedded vocoder: table cosine, fast exp2 / log2, 32 x 16 bit products (as dmr-vocoder's
   ambe_fx.h). */
#ifndef IMBE_FX_H
#define IMBE_FX_H

#include <stdint.h>

#include "imbe_tables.h"

/* One work area for the steps of a frame that never run at the same time (the synthesis' voiced and unvoiced parts) */
#define IMBE_SCRATCH_WORDS 1024
extern int32_t imbe_scratch[IMBE_SCRATCH_WORDS];

/* cos of a phase (2^32 = 2 pi), Q12: from the quarter table */
static inline int32_t cos_q12(uint32_t ph)
{
    uint32_t i = ph >> 20, q = i >> 10, j = i & 1023;
    int32_t v = ImbeCosQ12[(q & 1) ? 1024 - j : j];
    return ((q + 1) & 2) ? -v : v;
}

static inline int32_t iabs32(int32_t v)
{
    return v < 0 ? -v : v;
}

/* (a * b) >> 15 for a 32-bit a and a Q15 b, without a 64-bit product */
static inline int32_t mul32x16(int32_t a, int32_t b)
{
    return (a >> 15) * b + (((a & 0x7FFF) * b) >> 15);
}

/* log2 x, Q12, for x > 0 (x = 0 gives -32 << 12) */
static inline int32_t log2_q12(uint32_t x)
{
    if (x == 0) {
        return -(32 << 12);
    }
    int e = 31 - __builtin_clz(x);
    uint32_t m = x << (31 - e); /* 1.31 */
    uint32_t i = (m >> 23) & 0xFF, f = (m >> 15) & 0xFF;
    int32_t a = ImbeLog2Q12[i], b = ImbeLog2Q12[i + 1];
    return (e << 12) + a + (((b - a) * (int32_t)f) >> 8);
}

/* 2^(v / 4096) x 2^frac, saturated to 32 bits */
static inline uint32_t exp2_q12(int32_t v, int frac)
{
    int32_t i = v >> 12; /* floor */
    uint32_t f = (uint32_t)v & 0xFFF;
    uint32_t a = ImbeExp2Q15[f >> 4], b = ImbeExp2Q15[(f >> 4) + 1];
    uint32_t m = a + (((b - a) * (f & 15)) >> 4); /* 2^(f / 4096), Q15 */
    int sh = i + frac - 15;
    if (sh >= 16) {
        return 0xFFFFFFFFu;
    }
    if (sh >= 0) {
        return m << sh;
    }
    return sh <= -32 ? 0 : m >> -sh;
}

#endif
