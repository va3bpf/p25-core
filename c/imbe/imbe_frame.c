/*
 * IMBE frame layer: 144 channel bits <-> 88 parameter bits (TIA-102.BABA's error control: the Golay and Hamming
 * codes and the scrambling of c1..c6; the parameters' order as mbelib's mbe_eccImbe7200x4400Data).
 *
 * Golay(23,12): the cyclic code of g(x) = x^11 + x^10 + x^6 + x^5 + x^4 + x^2 + 1 (0xC75), systematic, the 12 data bits
 * high. Decoded by trying the error patterns of the data bits (a perfect code: every word is within 3 bits of exactly
 * one codeword). Hamming(15,11): the 4 parity checks below, data bits 14..4, parity 3..0; the bit to correct found
 * from the columns of the checks (mbelib's syndrome table maps some syndromes to the wrong bit).
 */
#include <string.h>

#include "imbe.h"

#define GOLAY_POLY 0xC75u

static int popcount(uint32_t x)
{
    return __builtin_popcount(x);
}

/* the 11 parity bits of 12 data bits: d(x) x^11 mod g(x) */
static uint32_t golay_parity(uint32_t d12)
{
    uint32_t r = d12 << 11;
    for (int b = 22; b >= 11; b--) {
        if (r & (1u << b)) {
            r ^= GOLAY_POLY << (b - 11);
        }
    }
    return r & 0x7FF;
}

uint32_t imbe_golay23_encode(uint32_t d12)
{
    d12 &= 0xFFF;
    return (d12 << 11) | golay_parity(d12);
}

uint32_t imbe_golay23_decode(uint32_t c23, int *errors)
{
    uint32_t d = (c23 >> 11) & 0xFFF, s = golay_parity(d) ^ (c23 & 0x7FF);
    uint32_t p[12];
    if (popcount(s) <= 3) {
        *errors = popcount(s);
        return d;
    }
    for (int i = 0; i < 12; i++) {
        p[i] = golay_parity(1u << i);
    }
    for (int i = 0; i < 12; i++) {
        uint32_t si = s ^ p[i];
        if (popcount(si) <= 2) {
            *errors = 1 + popcount(si);
            return d ^ (1u << i);
        }
    }
    for (int i = 0; i < 12; i++) {
        for (int j = i + 1; j < 12; j++) {
            uint32_t sij = s ^ p[i] ^ p[j];
            if (popcount(sij) <= 1) {
                *errors = 2 + popcount(sij);
                return d ^ (1u << i) ^ (1u << j);
            }
        }
    }
    for (int i = 0; i < 12; i++) {
        for (int j = i + 1; j < 12; j++) {
            for (int k = j + 1; k < 12; k++) {
                if ((s ^ p[i] ^ p[j] ^ p[k]) == 0) {
                    *errors = 3;
                    return d ^ (1u << i) ^ (1u << j) ^ (1u << k);
                }
            }
        }
    }
    *errors = -1; /* unreachable for a perfect code */
    return d;
}

/* the parity checks of Hamming(15,11), bit 14 first (row 0 gives syndrome bit 3) */
static const uint16_t HAM_H[4] = {0x7F08, 0x78E4, 0x66D2, 0x55B1};

static uint32_t hamming_syndrome(uint32_t c15)
{
    uint32_t s = 0;
    for (int r = 0; r < 4; r++) {
        s = (s << 1) | (uint32_t)(popcount(c15 & HAM_H[r]) & 1);
    }
    return s;
}

uint32_t imbe_hamming15_encode(uint32_t d11)
{
    uint32_t c = (d11 & 0x7FF) << 4;
    return c | hamming_syndrome(c); /* the parity bits are the syndrome of the data alone: each check has one */
}

uint32_t imbe_hamming15_decode(uint32_t c15, int *errors)
{
    uint32_t s = hamming_syndrome(c15);
    *errors = 0;
    if (s) {
        for (int p = 0; p < 15; p++) {
            if (hamming_syndrome(1u << p) == s) {
                c15 ^= 1u << p;
                *errors = 1;
                break;
            }
        }
    }
    return (c15 >> 4) & 0x7FF;
}

/* ---- bits ---- */

static uint32_t get_bits(const uint8_t *buf, int pos, int n)
{
    uint32_t v = 0;
    for (int i = 0; i < n; i++, pos++) {
        v = (v << 1) | ((buf[pos >> 3] >> (7 - (pos & 7))) & 1u);
    }
    return v;
}

static void put_bits(uint8_t *buf, int pos, int n, uint32_t v)
{
    for (int i = n - 1; i >= 0; i--, pos++) {
        uint8_t m = (uint8_t)(0x80 >> (pos & 7));
        if ((v >> i) & 1) {
            buf[pos >> 3] |= m;
        } else {
            buf[pos >> 3] &= (uint8_t)~m;
        }
    }
}

static uint32_t u_get(const uint8_t *u, int pos, int n)
{
    uint32_t v = 0;
    for (int i = 0; i < n; i++) {
        v = (v << 1) | (u[pos + i] & 1u);
    }
    return v;
}

static void u_put(uint8_t *u, int pos, int n, uint32_t v)
{
    for (int i = 0; i < n; i++) {
        u[pos + i] = (uint8_t)((v >> (n - 1 - i)) & 1);
    }
}

static const uint8_t CW_LEN[8] = {23, 23, 23, 23, 15, 15, 15, 7}; /* c0..c7 */
static const uint8_t U_LEN[8] = {12, 12, 12, 12, 11, 11, 11, 7};  /* u0..u7 */

/* the scrambling of c1..c6 (114 bits, c1 bit 22 first): p(0) = 16 u0, p(n) = 173 p(n-1) + 13849 mod 65536, the
   n-th bit (n = 1..114) the top bit of p(n). Returns the mask of word i (1..6), MSB first. */
static void pn_masks(uint32_t u0, uint32_t mask[8])
{
    uint32_t p = 16 * u0;
    for (int i = 1; i <= 6; i++) {
        uint32_t m = 0;
        for (int b = 0; b < CW_LEN[i]; b++) {
            p = (173 * p + 13849) & 0xFFFF;
            m = (m << 1) | (p >> 15);
        }
        mask[i] = m;
    }
}

void imbe_decode_frame(const uint8_t in18[18], uint8_t u88[88], imbe_errors *err)
{
    uint32_t c[8], u[8], mask[8];
    int pos = 0, e;
    for (int i = 0; i < 8; i++) {
        c[i] = get_bits(in18, pos, CW_LEN[i]);
        pos += CW_LEN[i];
    }
    err->e0 = err->etotal = err->bad = 0;
    u[0] = imbe_golay23_decode(c[0], &e);
    err->e0 = e;
    err->etotal = e;
    pn_masks(u[0], mask);
    for (int i = 1; i <= 6; i++) {
        uint32_t w = c[i] ^ mask[i];
        u[i] = i <= 3 ? imbe_golay23_decode(w, &e) : imbe_hamming15_decode(w, &e);
        err->etotal += e;
    }
    u[7] = c[7];
    pos = 0;
    for (int i = 0; i < 8; i++) {
        u_put(u88, pos, U_LEN[i], u[i]);
        pos += U_LEN[i];
    }
}

void imbe_encode_frame(const uint8_t u88[88], uint8_t out18[18])
{
    uint32_t u[8], c[8], mask[8];
    int pos = 0;
    for (int i = 0; i < 8; i++) {
        u[i] = u_get(u88, pos, U_LEN[i]);
        pos += U_LEN[i];
    }
    pn_masks(u[0], mask);
    c[0] = imbe_golay23_encode(u[0]);
    for (int i = 1; i <= 6; i++) {
        c[i] = (i <= 3 ? imbe_golay23_encode(u[i]) : imbe_hamming15_encode(u[i])) ^ mask[i];
    }
    c[7] = u[7];
    memset(out18, 0, 18);
    pos = 0;
    for (int i = 0; i < 8; i++) {
        put_bits(out18, pos, CW_LEN[i], c[i]);
        pos += CW_LEN[i];
    }
}
