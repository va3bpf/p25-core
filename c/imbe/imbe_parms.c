/*
 * IMBE parameter decoding, 88 bits -> w0, L, voicing and the log2 spectral amplitudes. Ported from mbelib's
 * mbe_decodeImbe4400Parms (ISC licence, see LICENSE), in fixed point: the gains, the DCT coefficients and log2 Ml in
 * Q12, Ml linear in Q6.
 *
 * Two differences from mbelib, both as the standard has them: the prediction's interpolation (eq. 75, 76) in exact
 * integers (mbelib's float can land on the wrong side of a whole number), and log2 M0(-1) = log2 M1(-1).
 */
#include <string.h>

#include "imbe.h"
#include "imbe_fx.h"

#define PH24 178956971u /* 2^32 / 24 */

void imbe_parms_init(imbe_parms *p)
{
    memset(p, 0, sizeof(*p));
    p->w0 = 0.09378f;
    p->w0q = 64107620; /* 0.09378 / 2 pi x 2^32 */
    p->L = 30;
    for (int l = 0; l < 58; l++) {
        p->PSIl[l] = 1u << 30; /* pi / 2 */
    }
}

int imbe_decode_parms(const uint8_t u[88], imbe_parms *cur, imbe_parms *prev)
{
    cur->repeat = prev->repeat;

    /* b0: u0's top 6 bits and u7's bits 2, 1 (imbe_d 85, 86) */
    int b0 = (u[0] << 7) | (u[1] << 6) | (u[2] << 5) | (u[3] << 4) | (u[4] << 3) | (u[5] << 2) | (u[85] << 1) | u[86];
    if (b0 > 207) {
        return 1;
    }
    int L = ImbeL[b0], L9 = L - 9;
    cur->w0q = ImbeW0q[b0];
    cur->w0 = (float)cur->w0q * 1.4629181e-9f; /* 2 pi / 2^32 */
    cur->L = L;
    int K = L < 37 ? (L + 2) / 3 : 12;

    /* the other 79 bits into b1 .. b(L+1) (bit order of eq. 66 .. 70) */
    uint32_t b[58];
    memset(b, 0, sizeof(b));
    for (int i = 6; i < 85; i++) {
        const uint8_t *o = ImbeBo[L9][i - 6];
        b[o[0]] |= (uint32_t)u[i] << o[1];
    }

    /* voicing: one bit a band of 3 harmonics, b1 MSB first, the last band takes the rest */
    for (int l = 1; l <= L; l++) {
        int band = (l - 1) / 3;
        if (band > K - 1) {
            band = K - 1;
        }
        cur->Vl[l] = (int)((b[1] >> (K - 1 - band)) & 1);
    }

    /* the gain vector: G1 from b2, G2..G6 from b3..b7 (eq. 68) */
    int32_t Gm[7];
    Gm[1] = ImbeB2Q12[b[2] & 63];
    for (int i = 2; i <= 6; i++) {
        int bits = ImbeBaBits[L9][i - 2];
        int32_t x = 2 * (int32_t)b[i + 1] - (1 << bits) + 1; /* 2 (bm - 2^(B-1) + 0.5) */
        Gm[i] = (int32_t)(((int64_t)x * ImbeBaStepQ24[L9][i - 2]) >> 13);
    }
    /* inverse DCT of the gains: Ri = C(i,1) */
    int32_t Cik[7][11];
    for (int i = 1; i <= 6; i++) {
        int32_t s = 0;
        for (int m = 1; m <= 6; m++) {
            s += ((m == 1 ? 1 : 2) * Gm[m] * cos_q12((uint32_t)((m - 1) * (2 * i - 1)) * PH24)) >> 12;
        }
        Cik[i][1] = s;
    }
    /* the higher order coefficients from b8 .. b(L+1) */
    int m = 8;
    for (int i = 1; i <= 6; i++) {
        for (int k = 2; k <= ImbeBlockLen[L9][i - 1]; k++, m++) {
            int bits = ImbeHoba[L9][m - 8];
            if (bits == 0) {
                Cik[i][k] = 0;
            } else {
                int32_t x = 2 * (int32_t)b[m] - (1 << bits) + 1;
                Cik[i][k] = (x * ImbeQsSdQ16[bits - 1][k - 2]) >> 5; /* x / 2, Q16 -> Q12 */
            }
        }
    }
    /* inverse DCT of each block: the prediction residuals Tl */
    int32_t Tl[57];
    int n = 1;
    for (int i = 1; i <= 6; i++) {
        int ji = ImbeBlockLen[L9][i - 1];
        uint32_t step = (1u << 30) / (uint32_t)ji; /* (k-1)(j-1/2)/(2 ji) turns = (k-1)(2j-1) step */
        for (int j = 1; j <= ji; j++) {
            int32_t s = 0;
            for (int k = 1; k <= ji; k++) {
                s += ((k == 1 ? 1 : 2) * Cik[i][k] * cos_q12((uint32_t)((k - 1) * (2 * j - 1)) * step)) >> 12;
            }
            Tl[n++] = s;
        }
    }

    /* log2 Ml: Tl plus rho times the previous frame's, interpolated to the new L, less its mean (eq. 75 .. 79) */
    int32_t rho = L <= 15 ? 1638 : L <= 24 ? ((3 * L - 5) * 4096 + 50) / 100 : 2867; /* 0.4, 0.03 L - 0.05, 0.7 */
    if (cur->L > prev->L) {
        for (int l = prev->L + 1; l <= cur->L; l++) {
            prev->Ml[l] = prev->Ml[prev->L];
            prev->log2Ml[l] = prev->log2Ml[prev->L];
        }
    }
    prev->log2Ml[0] = prev->log2Ml[1];
    prev->Ml[0] = prev->Ml[1];
    int32_t pred[57], psum = 0;
    for (int l = 1; l <= L; l++) {
        int t = prev->L * l;
        int k = t / L;
        int32_t frac = ((t % L) << 12) / L; /* delta, Q12 */
        int32_t a = prev->log2Ml[k], c = prev->log2Ml[k + 1];
        pred[l] = a + (((c - a) * frac) >> 12);
        psum += pred[l];
    }
    int32_t sum77 = (int32_t)(((int64_t)rho * psum / L) >> 12);
    for (int l = 1; l <= L; l++) {
        int32_t v = Tl[l] + ((rho * pred[l]) >> 12) - sum77;
        cur->log2Ml[l] = v;
        uint32_t ml = exp2_q12(v, 6);
        cur->Ml[l] = (int32_t)(ml > 0xFFFFFF ? 0xFFFFFF : ml); /* 2^18: far above 16-bit PCM */
    }
    return 0;
}
