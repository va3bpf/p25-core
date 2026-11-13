/*
 * IMBE encoder: the parameter quantizer (TIA-102.BAAA's 7200x4400 frame, the inverse of imbe_parms.c), fixed point.
 *
 * In: the fundamental as b0, a voicing decision a harmonic and log2 of the spectral amplitudes (Q12). The quantizer
 * keeps a mirror of the decoder's state (its last decoded frame: the prediction of the amplitudes runs on what the
 * decoder has, not on what was analysed), and every decision is the nearest of the decoder's reconstruction values,
 * so a frame decoded and quantized again gives the same bits.
 */
#include <string.h>

#include "imbe.h"
#include "imbe_fx.h"

#define PH24 178956971u /* 2^32 / 24: (m - 1)(2 i - 1) / 24 turns */

void imbe_quant_init(imbe_quantizer *q)
{
    imbe_parms_init(&q->prev);
}

int imbe_b0_of(uint32_t w0q)
{
    /* ImbeW0q falls as b0 rises: the nearest */
    int lo = 0, hi = 207;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (ImbeW0q[mid] > w0q) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo > 0 && (ImbeW0q[lo - 1] - w0q) < (w0q - ImbeW0q[lo])) {
        lo--;
    }
    return lo;
}

static int clampi(int v, int lo, int hi)
{
    return v < lo ? lo : v > hi ? hi : v;
}

/* the nearest b of a uniform quantizer whose values are (b - 2^(B-1) + 1/2) step: b = floor(x / step + 2^(B-1)), with
   num / den = x / step */
static int uniform(int64_t num, int64_t den, int bits)
{
    int64_t half = (int64_t)1 << (bits - 1);
    int64_t v = num >= 0 ? num / den : -((-num + den - 1) / den); /* floor */
    return clampi((int)(v + half), 0, (1 << bits) - 1);
}

void imbe_quantize(imbe_quantizer *q, int b0, const uint8_t *Vl, const int32_t *log2Ml, uint8_t u[88])
{
    imbe_parms *prev = &q->prev;
    int L = ImbeL[b0], L9 = L - 9;
    int K = L < 37 ? (L + 2) / 3 : 12;
    uint32_t b[58];
    memset(b, 0, sizeof(b));
    b[0] = (uint32_t)b0;

    /* b1: a band of 3 harmonics voiced when most of its harmonics are (the last band takes the rest) */
    for (int k = 0; k < K; k++) {
        int first = 3 * k + 1, last = (k == K - 1) ? L : 3 * k + 3, v = 0;
        for (int l = first; l <= last; l++) {
            v += Vl[l] ? 1 : 0;
        }
        if (2 * v >= last - first + 1) {
            b[1] |= 1u << (K - 1 - k);
        }
    }

    /* the prediction from the decoder's last frame (as imbe_decode_parms, on a copy: it extends prev's amplitudes) */
    int32_t pl[58];
    memcpy(pl, prev->log2Ml, sizeof(pl));
    for (int l = prev->L + 1; l <= L; l++) {
        pl[l] = pl[prev->L];
    }
    pl[0] = pl[1];
    int32_t rho = L <= 15 ? 1638 : L <= 24 ? ((3 * L - 5) * 4096 + 50) / 100 : 2867;
    int32_t pred[57], psum = 0;
    for (int l = 1; l <= L; l++) {
        int t = prev->L * l;
        int k = t / L;
        int32_t frac = ((t % L) << 12) / L;
        int32_t a = pl[k], c = pl[k + 1];
        pred[l] = a + (((c - a) * frac) >> 12);
        psum += pred[l];
    }
    int32_t sum77 = (int32_t)(((int64_t)rho * psum / L) >> 12);
    /* the residuals the decoder must rebuild: log2 Ml = Tl + rho pred - sum77 */
    int32_t Tl[57];
    for (int l = 1; l <= L; l++) {
        Tl[l] = log2Ml[l] - ((rho * pred[l]) >> 12) + sum77;
    }

    /* each block's DCT: C(i,k) = 1/ji sum_j Tl cos(pi (k-1)(j-1/2) / ji) */
    int32_t Cik[7][11];
    int n = 1;
    for (int i = 1; i <= 6; i++) {
        int ji = ImbeBlockLen[L9][i - 1];
        uint32_t step = (1u << 30) / (uint32_t)ji;
        for (int k = 1; k <= ji; k++) {
            int64_t s = 0;
            for (int j = 1; j <= ji; j++) {
                s += (int64_t)Tl[n + j - 1] * cos_q12((uint32_t)((k - 1) * (2 * j - 1)) * step);
            }
            Cik[i][k] = (int32_t)((s / ji) >> 12);
        }
        n += ji;
    }
    /* the gain vector: G(m) = 1/6 sum_i C(i,1) cos(pi (m-1)(i-1/2) / 6) */
    int32_t Gm[7];
    for (int m = 1; m <= 6; m++) {
        int64_t s = 0;
        for (int i = 1; i <= 6; i++) {
            s += (int64_t)Cik[i][1] * cos_q12((uint32_t)((m - 1) * (2 * i - 1)) * PH24);
        }
        Gm[m] = (int32_t)((s / 6) >> 12);
    }
    /* G1: the nearest of b2's 64 values */
    int best = 0;
    for (int v = 1; v < 64; v++) {
        if (iabs32(ImbeB2Q12[v] - Gm[1]) < iabs32(ImbeB2Q12[best] - Gm[1])) {
            best = v;
        }
    }
    b[2] = (uint32_t)best;
    /* G2..G6: uniform, Gm = (b - 2^(B-1) + 1/2) step, step = StepQ24 / 2^12 in Q12 */
    for (int i = 2; i <= 6; i++) {
        int bits = ImbeBaBits[L9][i - 2];
        b[i + 1] = (uint32_t)uniform((int64_t)Gm[i] << 12, ImbeBaStepQ24[L9][i - 2], bits);
    }
    /* the higher order coefficients: C = (b - 2^(B-1) + 1/2) QsSd, QsSd = QsSdQ16 / 16 in Q12 */
    int m = 8;
    for (int i = 1; i <= 6; i++) {
        for (int k = 2; k <= ImbeBlockLen[L9][i - 1]; k++, m++) {
            int bits = ImbeHoba[L9][m - 8];
            if (bits) {
                b[m] = (uint32_t)uniform((int64_t)Cik[i][k] << 4, ImbeQsSdQ16[bits - 1][k - 2], bits);
            }
        }
    }

    /* into the 88 bits: b0 in u0's top 6 and u7's bits 2, 1; the rest by the bit order */
    memset(u, 0, 88);
    for (int i = 0; i < 6; i++) {
        u[i] = (uint8_t)((b0 >> (7 - i)) & 1);
    }
    u[85] = (uint8_t)((b0 >> 1) & 1);
    u[86] = (uint8_t)(b0 & 1);
    for (int i = 6; i < 85; i++) {
        const uint8_t *o = ImbeBo[L9][i - 6];
        u[i] = (uint8_t)((b[o[0]] >> o[1]) & 1);
    }

    /* the decoder's view of this frame, for the next prediction */
    imbe_parms cur = *prev;
    imbe_decode_parms(u, &cur, prev);
    *prev = cur;
}
