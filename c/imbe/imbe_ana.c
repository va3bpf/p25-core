/*
 * IMBE encoder: 160 samples -> the 88 parameter bits (imbe_enc_frame), or 18 bytes with the error-correction
 * (imbe_enc_frame18).
 *
 * The analysis is dmr-vocoder's AMBE+2 encoder's (c/ambe/ambe_enc.c, the same author; IMBE and AMBE+2 are both MBE
 * vocoders and analyse speech the same way), with IMBE's parameters: its pitch grid (b0 0..207), voicing a band of
 * three harmonics, no silence frame (quiet frames are unvoiced noise); then the quantizer (imbe_enc.c).
 *
 * Frame k is analysed around input sample 160 k - 200. The pitch tracker looks two frames ahead (E(P) of frame
 * k + 2 reaches sample 160 k + 280), so frame k is encoded once frame k + 1 is in: one frame of delay.
 *
 * Integer, for the CK803S: the spectrum is a fixed-point FFT scaled to 12 bits, the harmonic fit sums in 32 bits,
 * E(P) and the tracker in Q16, amplitudes in log2 Q12. Float only where a frame needs a few dozen operations (the
 * refinement's grid, the voice activity in dB).
 */
#include <string.h>

#include "imbe.h"
#include "imbe_ana_fx.h"
#include "imbe_tables.h"

#define PROF(i, x) x

int32_t imbea_scratch[IMBEA_SCRATCH_WORDS];

#define N IMBE_N
#define NW 256
#ifndef IMBEA_FSH
#define IMBEA_FSH 1               /* the spectrum's FFT: 1024 points >> IMBEA_FSH (512: as good, half the work) */
#endif
#define FSH IMBEA_FSH
#define NFFT (1024 >> FSH)
#define NBINS (NFFT / 2 + 1)
#define TSH (11 - FSH)           /* the window transform table has 2^(16 - TSH) points a bin */
#define CENTER (-200)
#define NI 301
#define NP 217                   /* the pitch grid: P = 20 .. 128 in half samples */
#define TRACK_PMAX_I 160         /* grid index of P = 100: longer periods count as aperiodic */
#define Q16(v) ((int32_t)((v) * 65536.0 + ((v) < 0 ? -0.5 : 0.5)))
#define Q12(v) ((int32_t)((v) * 4096.0 + ((v) < 0 ? -0.5 : 0.5)))
/* the constants of py/ambe_analysis.py */
#define VOICED_ERR 0.8f
#define VOICING_HYSTERESIS 2.0f
#define REFINE_FMAX 3000.0f
#define QUIET_FULL (-25.0f)
#define QUIET_NONE (-42.0f)
#define HOLD_E Q16(0.8)
#define VOICE_MAX_E 0.85f
#define VOICE_E_FULL 0.4f
#define VOICE_E_NONE 0.8f
#define LONG_P 95.0f
#define LONG_MAX_E 0.3f
#define SUBMULT_TOL Q16(0.05)
#define SUBMULT_MAXERR Q16(1.0)
/* py/imbe_encoder.py: voice activity */
#define VAD_MARGIN 15.0f
#define VAD_MIN_DB 20.0f
#define FLOOR_RISE 0.05f
#define PEAK_DECAY 0.05f
#define VAD_BELOW_PEAK 45.0f
/* log2 constants, Q12: amplitudes from the fit (see harmonics()) */
/* AMBE's (its decoder's output gain 1.7 included) moved to IMBE's level by these, Q12: measured against the
   reference encoder's amplitudes on speech (py/enc_test.py) */
#ifndef IMBE_ENC_V_OFS
#define IMBE_ENC_V_OFS 0 /* the voiced as AMBE's */
#endif
#ifndef IMBE_ENC_U_OFS
#define IMBE_ENC_U_OFS (-6300) /* -1.54: the unvoiced 9 dB under (PESQ the same from 0 / -4300 to here; the level 0.7 dB over the input) */
#endif
#define LOG2_AMP_V (21005 + IMBE_ENC_V_OFS)         /* 12 - log2(sum of the window) + 1 - log2(OUT_GAIN 1.7) */
#define LOG2_AMP_U (-33166 + (FSH << 11) + IMBE_ENC_U_OFS) /* 1 - log2(NFFT sum(window^2)) / 2 - log2(1.7) */
#define LOG2_FLOOR (-40820)      /* log2(1e-3) */

static int bits(uint32_t v) /* the number of bits of v: 0 for 0 */
{
    return v ? 32 - __builtin_clz(v) : 0;
}

static inline int32_t iabs(int32_t v)
{
    return v < 0 ? -v : v;
}

/* (a * b) >> 12, |a| < 2^28, |b| < 2^15, no 64-bit product */
static inline int32_t mulq12s(int32_t a, int32_t b)
{
    return (a >> 12) * b + (((a & 0xFFF) * b) >> 12);
}

static inline int32_t mulq14s(int32_t a, int32_t b)
{
    return (a >> 14) * b + (((a & 0x3FFF) * b) >> 14);
}

static inline int32_t mulq15s(int32_t a, int32_t b)
{
    return (a >> 15) * b + (((a & 0x7FFF) * b) >> 15);
}

static int32_t log2_u64(uint64_t v) /* Q12 */
{
    uint32_t hi = (uint32_t)(v >> 32);
    if (hi) {
        int s = bits(hi);
        return log2_q12((uint32_t)(v >> s)) + (s << 12);
    }
    return log2_q12((uint32_t)v);
}

/* ---- input ---- */

static uint32_t spec_e[NBINS + 1]; /* sum of |S|^2 over the bins below k (mod 2^32: a band's sum fits 32 bits) */

static int32_t xin(const imbe_encoder *e, long a)
{
    long i = a - (e->n_in - IMBE_ENC_HIST);
    return i < 0 || i >= IMBE_ENC_HIST ? 0 : e->x[i];
}

/* ---- spectrum: the windowed frame (mean removed), a 1024-point FFT, the window's linear phase taken out, scaled
   to 12 bits: S = sr + j si, true value S x 2^ex ---- */

/* in-place NFFT / 2-point complex FFT (e^-j), Q15 twiddles; the data stay within 15 bits (in: |x| < 2^14), every
   stage halves them: out = FFT / (NFFT / 2) */
static IMBE_RAM void fft_half(int32_t *re, int32_t *im)
{
    const int n = NFFT / 2;
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            int32_t t = re[i];
            re[i] = re[j];
            re[j] = t;
            t = im[i];
            im[i] = im[j];
            im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int step = 2048 / len; /* e^(-j 2 pi k / len) = cos(2 pi k step / 2048) - j sin */
        for (int k = 0; k < len / 2; k++) {
            int32_t wr = cos2048_q15(k * step), wi = -cos2048_q15(k * step - 512);
            for (int i = k; i < n; i += len) {
                int32_t *ar = &re[i], *ai = &im[i], *br = &re[i + len / 2], *bi = &im[i + len / 2];
                int32_t vr = (*br * wr - *bi * wi) >> 15, vi = (*br * wi + *bi * wr) >> 15;
                *br = (*ar - vr) >> 1;
                *bi = (*ai - vi) >> 1;
                *ar = (*ar + vr) >> 1;
                *ai = (*ai + vi) >> 1;
            }
        }
    }
}

/* the frame's spectrum: the 1024-point FFT of NW real samples as a 512-point complex one of the even and odd
   samples, X[k] = (Z[k] + Z*[512-k]) / 2 + e^(-j 2 pi k / 1024) (Z[k] - Z*[512-k]) / 2j; the window's linear phase
   taken out (x e^(j w (NW - 1) / 2)), scaled to 12 bits. S = sr + j si, true value S x 2^ex */
static int spectrum(const imbe_encoder *e, long k, int16_t *sr, int16_t *si)
{
    int32_t *re = imbea_scratch, *im = imbea_scratch + NFFT / 2, *xr = im + NFFT / 2, *xi = xr + NBINS;
    long a = N * k + CENTER - NW / 2;
    int32_t mean = 0, mx = 0, v[NW];
    for (int i = 0; i < NW; i++) {
        v[i] = xin(e, a + i);
        mean += v[i];
    }
    mean /= NW;
    for (int i = 0; i < NW; i++) {
        v[i] = ((v[i] - mean) * ImbeAHammingQ15[i]) >> 15;
        mx = iabs(v[i]) > mx ? iabs(v[i]) : mx;
    }
    int g = 13 - bits((uint32_t)mx); /* the frame to 13 bits (the FFT halves it 9 - FSH times) */
    memset(re, 0, sizeof(int32_t) * NFFT);
    for (int i = 0; i < NW / 2; i++) {
        re[i] = g >= 0 ? v[2 * i] << g : v[2 * i] >> -g;
        im[i] = g >= 0 ? v[2 * i + 1] << g : v[2 * i + 1] >> -g;
    }
    PROF(8, fft_half(re, im));
    mx = 0;
    for (int i = 0; i < NBINS; i++) {
        int j = (NFFT / 2 - i) & (NFFT / 2 - 1), ii = i & (NFFT / 2 - 1);
        int32_t ar = re[ii], ai = im[ii], br = re[j], bi = -im[j];
        int32_t er = (ar + br) >> 1, ei = (ai + bi) >> 1, odr = (ai - bi) >> 1, odi = -((ar - br) >> 1);
        int tw = (2 * i) << FSH; /* e^(-j 2 pi i / NFFT) */
        int32_t wr = cos2048_q15(tw), wi = -cos2048_q15(tw - 512);
        int32_t x_r = er + ((odr * wr - odi * wi) >> 15), x_i = ei + ((odr * wi + odi * wr) >> 15);
        int idx = ((255 * i) << FSH) & 2047; /* x e^(j 2 pi i 127.5 / NFFT) */
        int32_t c = cos2048_q15(idx), sn = cos2048_q15(idx - 512);
        xr[i] = (x_r * c - x_i * sn) >> 15;
        xi[i] = (x_r * sn + x_i * c) >> 15;
        mx = iabs(xr[i]) > mx ? iabs(xr[i]) : mx;
        mx = iabs(xi[i]) > mx ? iabs(xi[i]) : mx;
    }
    int h = bits((uint32_t)mx) - 12; /* to 12 bits */
    for (int i = 0; i < NBINS; i++) {
        sr[i] = (int16_t)(h >= 0 ? xr[i] >> h : xr[i] << -h);
        si[i] = (int16_t)(h >= 0 ? xi[i] >> h : xi[i] << -h);
        spec_e[i + 1] = spec_e[i] + (uint32_t)(sr[i] * sr[i] + si[i] * si[i]);
    }
    return h - g + 9 - FSH;
}

static float frame_energy_db(const imbe_encoder *e, long k)
{
    long a = N * k + CENTER - NW / 2;
    int32_t mean = 0;
    for (int i = 0; i < NW; i++) {
        mean += xin(e, a + i);
    }
    mean /= NW;
    uint64_t p = 0;
    for (int i = 0; i < NW; i++) {
        int32_t v = ((xin(e, a + i) - mean) * ImbeAHammingQ15[i]) >> 15;
        p += (uint32_t)(v * v);
    }
    if (p == 0) {
        return -30.0f; /* 10 log10(1e-3) */
    }
    /* 10 log10(p / sum(window^2)); log2(101.3434) = 6.66310 */
    return (log2_u64(p) - Q12(6.66310)) * (3.0103f / 4096);
}

/* ---- the harmonic fit ----
   For each harmonic l of w0 (wb: w0 in bins, Q16) over the bins (l - 1/2) w0 <= w_k < (l + 1/2) w0: the window's
   transform W at the bin's distance from l w0 (real: spectrum() took out its phase), nr + j ni = sum W S,
   ww = sum W^2, en = sum |S|^2 and proj = |sum W S|^2 / ww, the energy the fitted sinusoid explains. */

typedef struct {
    uint32_t en[58], proj[58], ww[58];
} fit_t;

#define WINT ImbeAWinTQ12

static fit_t fit; /* the fit of the moment: each user is done with it before the next fit */

static uint32_t projection(int32_t nr, int32_t ni, uint32_t ww)
{
    if (ww == 0) {
        return 0;
    }
    int32_t ax = iabs(nr) > iabs(ni) ? iabs(nr) : iabs(ni);
    if (ax == 0) {
        return 0;
    }
    int sn = bits((uint32_t)ax) - 15; /* nr, ni to 15 bits */
    int32_t a = sn >= 0 ? nr >> sn : nr << -sn, b = sn >= 0 ? ni >> sn : ni << -sn;
    uint32_t num = (uint32_t)(a * a) + (uint32_t)(b * b);
    int t = bits(ww) - 16;
    t = t < 0 ? 0 : t;
    uint32_t q = num / (ww >> t);
    int ex = 2 * sn - t;
    if (ex >= 0) {
        return ex >= 31 || q > (0xFFFFFFFFu >> ex) ? 0xFFFFFFFFu : q << ex;
    }
    return ex <= -32 ? 0 : q >> -ex;
}

static IMBE_RAM void harmonic_fit(const int16_t *sr, const int16_t *si, uint32_t wb, int L, fit_t *f)
{
#ifdef IMBEA_PROFILE
    imbea_prof[16]++;         /* fits */
    imbea_prof[17] += L;      /* harmonics fitted */
#endif
    for (int l = 1; l <= L; l++) {
        uint32_t c = wb * (uint32_t)l, lo = c - wb / 2, hi = c + wb / 2;
        int k0 = (int)((lo + 0xFFFF) >> 16), k1 = (int)((hi - 1) >> 16);
        int32_t nr = 0, ni = 0;
        uint32_t ww = 0;
        if (k1 > NBINS - 1) {
            k1 = NBINS - 1;
        }
        uint32_t en = k1 >= k0 ? spec_e[k1 + 1] - spec_e[k0] : 0;
        /* W at the nearest of 64 points a bin from l w0; a bin further is 64 points further. The bins of a
           harmonic are within w0 / 2 <= 27 bins: never past the table. */
        int kc = (int)(c >> 16), k = k0;
        const short *w = WINT + ((c - ((uint32_t)k << 16) + (1u << (TSH - 1))) >> TSH);
        for (; k <= k1 && k <= kc; k++, w -= 1 << (16 - TSH)) { /* below l w0 */
            nr += *w * sr[k];
            ni += *w * si[k];
            ww += (uint32_t)(*w * *w);
        }
        if (k <= k1) {
            w = WINT + ((((uint32_t)k << 16) - c + (1u << (TSH - 1))) >> TSH);
            for (; k <= k1; k++, w += 1 << (16 - TSH)) { /* above */
                nr += *w * sr[k];
                ni += *w * si[k];
                ww += (uint32_t)(*w * *w);
            }
        }
        uint32_t p = projection(nr, ni, ww);
        f->en[l] = en;
        f->ww[l] = ww;
        f->proj[l] = p < en ? p : en;
    }
}

static uint32_t bins_of_period(float q) /* w0 = 2 pi / q in bins, Q16 */
{
    return (uint32_t)((67108864 >> FSH) / q);
}

/* sum(err) / sum(en) below REFINE_FMAX, Q16 */
static int32_t fit_error(const int16_t *sr, const int16_t *si, float q)
{
    fit_t *f = &fit;
    int L = (int)(REFINE_FMAX / (8000.0f / q));
    L = L < 1 ? 1 : L;
    PROF(5, harmonic_fit(sr, si, bins_of_period(q), L, f));
    uint64_t se = 0, sn = 0;
    for (int l = 1; l <= L; l++) {
        se += f->en[l] - f->proj[l];
        sn += f->en[l];
    }
    return sn ? (int32_t)((se << 16) / sn) : 0;
}

/* ---- pitch tracker (py Analyzer.E, track) ---- */

static IMBE_RAM void compute_E(const imbe_encoder *e, long k, int32_t *E)
{
    long a = N * k + CENTER - NI / 2;
    int16_t *xs = (int16_t *)imbea_scratch, *ss = xs + NI, *y = ss + NI; /* 2 x 301 int16, y + 160 zeros */
    int32_t *r = imbea_scratch + 540;                                     /* after them, 152 words */
    int16_t *xw = (int16_t *)(imbea_scratch + 700);                       /* the input, 321 samples */
    int32_t mx = 0;
    for (int i = 0; i < NI + 20; i++) {
        xw[i] = (int16_t)xin(e, a - 10 + i);
    }
    for (int i = 0; i < NI; i++) {
        int32_t s = ImbeALpfQ15[10] * xw[i + 10];
        for (int t = 0; t < 10; t++) { /* the filter is symmetric */
            s += ImbeALpfQ15[t] * (xw[i + t] + xw[i + 20 - t]);
        }
        ss[i] = (int16_t)(s >> 15);
        xs[i] = xw[i + 10];
        mx = iabs(xs[i]) > mx ? iabs(xs[i]) : mx;
        mx = iabs(ss[i]) > mx ? iabs(ss[i]) : mx;
    }
    int g = 11 - bits((uint32_t)mx); /* to 11 bits */
    int32_t n0 = 0, raw = 0;
    for (int i = 0; i < NI && mx; i++) {
        int32_t s = g >= 0 ? ss[i] << g : ss[i] >> -g, x = g >= 0 ? xs[i] << g : xs[i] >> -g;
        y[i] = (int16_t)((s * ImbeAWi2Q15[i]) >> 15);
        n0 += s * y[i];
        raw += x * ((x * ImbeAWi2Q15[i]) >> 15);
    }
    /* no periodicity: silence, or less than 3 % of the energy below 1 kHz (a fricative) */
    if (mx == 0 || n0 <= 0 || (int64_t)n0 * 100 < (int64_t)raw * 3) {
        for (int i = 0; i < NP; i++) {
            E[i] = 65536;
        }
        return;
    }
    memset(y + NI, 0, 160 * sizeof(int16_t));
    r[0] = 0;
    for (int i = 0; i < NI; i++) {
        r[0] += y[i] * y[i];
    }
    for (int t = 20; t <= NI / 2 + 1; t += 4) { /* the lags of P >= 20, 20 .. 151, four at a time */
        int32_t c0 = 0, c1 = 0, c2 = 0, c3 = 0;
        const int16_t *p = y + t;
        for (int i = 0; i + t < NI; i++) {
            int32_t v = y[i];
            c0 += v * p[i];
            c1 += v * p[i + 1];
            c2 += v * p[i + 2];
            c3 += v * p[i + 3];
        }
        r[t] = c0;
        r[t + 1] = c1;
        r[t + 2] = c2;
        r[t + 3] = c3;
    }
    int sh = bits((uint32_t)r[0]) - 26; /* tot below fits 31 bits */
    if (sh > 0) {
        for (int t = 0; t <= NI / 2 + 1; t++) {
            r[t] >>= sh;
        }
        n0 >>= sh;
    }
    int32_t d = n0 >> 12; /* a 32-bit division a period: the ratio to 1 / 4096 */
    d = d < 1 ? 1 : d;
    for (int i = 0; i < NP; i++) {
        if (i > TRACK_PMAX_I) {
            E[i] = 65536;
            continue;
        }
        int p2 = 40 + i, nl = 300 / p2; /* 2 P; int(150 / P) */
        int32_t tot = 0;
        for (int m = 1; m <= nl; m++) {
            int lag2 = p2 * m, t = lag2 >> 1;
            tot += (lag2 & 1) ? (r[t] + r[t + 1]) >> 1 : r[t];
        }
        tot = r[0] + 2 * tot;
        /* E = (1 - P max(w^2) tot / n0) / (1 - P sum(w^4)) */
        int32_t rho = (mulq14s(tot, ImbeAPmQ14[i]) / d) << 4; /* (P max(w^2) tot / n0) in Q16 */
        rho = rho > (1 << 18) ? (1 << 18) : rho < -(1 << 18) ? -(1 << 18) : rho;
        E[i] = (((65536 - rho) >> 2) * ImbeAInvDQ12[i]) >> 10;
    }
}

/* E(P) of frame f, kept in slot f % 3 (the tracker wants frames k .. k + 2) */
static const int32_t *E_of(imbe_encoder *e, long f)
{
    int slot = (int)(f % 3);
    if (e->E_tag[slot] != f) {
        PROF(6, compute_E(e, f, e->E[slot]));
        e->E_tag[slot] = f;
    }
    return e->E[slot];
}

/* min over the +-20 % around P_j of cost + |ln P - ln P_j| (PITCH_JUMP_COST 1), and where */
static IMBE_RAM int32_t best_from(const int32_t *cost, int j, int *idx)
{
    int32_t best = 0x7FFFFFFF;
    int bi = j;
    int32_t lj = ImbeALogPQ16[j];
    for (int i = ImbeAPLo[j]; i < ImbeAPHi[j]; i++) {
        int32_t r = cost[i] + iabs(ImbeALogPQ16[i] - lj);
        if (r < best) {
            best = r;
            bi = i;
        }
    }
    if (idx) {
        *idx = bi;
    }
    return best;
}

/* best_from() for every j at once: min over [lo_j, j] of (cost_i - x_i) + x_j and over [j, hi_j) of (cost_i + x_i) - x_j,
   x = ln P; both windows only move up with j, so a monotone queue gives each minimum (wk: 4 NP words) */
static IMBE_RAM void best_all(const int32_t *cost, int32_t *out, int32_t *wk)
{
    int32_t *a = wk, *b = wk + NP;
    int *qa = (int *)(wk + 2 * NP), *qb = qa + NP;
    int ha = 0, ta = 0, hb = 0, tb = 0, nxt = 0;
    for (int j = 0; j < NP; j++) {
        a[j] = cost[j] - ImbeALogPQ16[j];
        while (ta > ha && a[qa[ta - 1]] >= a[j]) {
            ta--;
        }
        qa[ta++] = j;
        while (qa[ha] < ImbeAPLo[j]) {
            ha++;
        }
        while (nxt < ImbeAPHi[j]) {
            b[nxt] = cost[nxt] + ImbeALogPQ16[nxt];
            while (tb > hb && b[qb[tb - 1]] >= b[nxt]) {
                tb--;
            }
            qb[tb++] = nxt++;
        }
        while (qb[hb] < j) {
            hb++;
        }
        int32_t l = a[qa[ha]] + ImbeALogPQ16[j], r = b[qb[hb]] - ImbeALogPQ16[j];
        out[j] = l < r ? l : r;
    }
}

/* c <= num / den x cF, as py's c / max(cF, 1e-9) <= num / den */
static int ratio_le(int32_t c, int32_t cF, int num, int den)
{
    int64_t d = cF > 0 ? cF : 1;
    return (int64_t)c * den <= d * num;
}

static float track(imbe_encoder *e, long k)
{
    const int32_t *E0 = E_of(e, k), *E1 = E_of(e, k + 1), *E2 = E_of(e, k + 2);
    int32_t *c12 = imbea_scratch, *ceF = imbea_scratch + NP, *wk = imbea_scratch + 2 * NP; /* compute_E is done */
    int i1 = e->hist_i[0];
    int32_t e1 = e->hist_e[0], e2 = e->hist_e[1];
    int iB;
    int32_t ceB = best_from(E0, i1, &iB) + e1 + e2;
    best_all(E2, c12, wk);
    for (int j = 0; j < NP; j++) {
        c12[j] += E1[j];
    }
    best_all(c12, ceF, wk);
    int iF = 0;
    for (int i = 0; i < NP; i++) {
        ceF[i] += E0[i];
        if (ceF[i] < ceF[iF]) {
            iF = i;
        }
    }
    /* sub-multiples of the look-ahead period that are nearly as good */
    int32_t cF = ceF[iF];
    for (int d = (40 + iF) / 40; d > 1; d--) { /* int(P / 20) */
        int j = ((40 + iF) * 2 / d + 1) / 2 - 40; /* the grid point nearest P / d */
        j = j < 0 ? 0 : j;
        int32_t c = ceF[j];
        if ((c <= Q16(0.85) && ratio_le(c, cF, 17, 10)) || (c <= Q16(0.4) && ratio_le(c, cF, 35, 10)) ||
            (c <= Q16(0.05) && ratio_le(c, cF, 5, 1))) {
            iF = j;
            break;
        }
    }
    int i = (ceB <= Q16(0.48) || ceB <= ceF[iF] + iabs(ImbeALogPQ16[iF] - ImbeALogPQ16[i1])) ? iB : iF;
    e->hist_i[1] = i1;
    e->hist_e[1] = e1;
    e->hist_i[0] = i;
    e->hist_e[0] = E0[i];
    e->track_e = E0[i];
    return 20.0f + 0.5f * i;
}

/* ---- pitch refinement (py Analyzer.refine, pitch) ---- */

static void tmin(int32_t *e, float *q, int32_t e2, float q2) /* tuple min of (e, q) */
{
    if (e2 < *e || (e2 == *e && q2 < *q)) {
        *e = e2;
        *q = q2;
    }
}

static float refine(const int16_t *sr, const int16_t *si, float p, int32_t *err)
{
    const float span = 0.03f;
    int32_t e = 0;
    float q = 0;
    for (int i = 0; i < 13; i++) {
        float qi = p * (1 - span) + (2 * p * span) * i / 12;
        int32_t ei = fit_error(sr, si, qi);
        if (i == 0) {
            e = ei;
            q = qi;
        } else {
            tmin(&e, &q, ei, qi);
        }
    }
    float step = 0.0025f * p;
    for (int it = 0; it < 4; it++) {
        float q0 = q;
        int32_t em = fit_error(sr, si, q0 - step), ep = fit_error(sr, si, q0 + step);
        tmin(&e, &q, em, q0 - step);
        tmin(&e, &q, ep, q0 + step);
        step /= 2;
    }
    int32_t el = fit_error(sr, si, q - step), er = fit_error(sr, si, q + step);
    int32_t den = el - 2 * e + er;
    if (den > 0) {
        int32_t t = ((el - er) << 7) / den; /* 0.5 (el - er) / den, Q8 */
        t = t > 256 ? 256 : t < -256 ? -256 : t;
        q += step * t * (1.0f / 256);
    }
    int32_t m = e < el ? e : el;
    *err = m < er ? m : er;
    return q;
}

/* -> w0 (rad / sample) */
static float pitch(imbe_encoder *e, long k, const int16_t *sr, const int16_t *si)
{
    float p;
    PROF(4, p = track(e, k));
    int32_t err;
    float q = refine(sr, si, p, &err);
    int periodic = err < SUBMULT_MAXERR;
    for (int m = 3; m >= 2 && periodic; m--) {
        float pm = q / m;
        if (pm >= 20) {
            int32_t ek;
            float qk = refine(sr, si, pm, &ek);
            if (ek < err + SUBMULT_TOL) {
                err = ek;
                q = qk;
                break;
            }
        }
    }
    if (2 * q <= 100) {
        /* the second harmonic taken for the pitch: twice the period when its odd harmonics have energy */
        fit_t *f = &fit;
        float h = 2 * q;
        int n = (int)(2000 / (8000 / h));
        harmonic_fit(sr, si, bins_of_period(h), n, f);
        int lmin = (int)(400 / (8000 / h)); /* harmonics from 400 Hz */
        if (lmin * (8000 / h) < 400) {
            lmin++;
        }
        uint64_t odd = 0, even = 0;
        for (int l = lmin < 1 ? 1 : lmin; l <= n; l++) {
            uint64_t p2 = f->ww[l] ? ((uint64_t)f->proj[l] << 16) / f->ww[l] : 0; /* |A_l|^2 */
            if (l & 1) {
                odd += p2;
            } else {
                even += p2;
            }
        }
        /* nearer the last frame's pitch: |ln(h / last)| < |ln(w0 / last)|, as periods */
        int near = 0;
        if (e->have_last) {
            float lp = 6.2831853f / e->last_w0, r1 = h / lp, r2 = q / lp;
            r1 = r1 < 1 ? 1 / r1 : r1;
            r2 = r2 < 1 ? 1 / r2 : r2;
            near = r1 < r2;
        }
        if (near && even > 0 && odd * 1000 > even * 251) { /* within 6 dB */
            q = h;
        }
    }
    float w0 = 6.2831853f / q;
    if (e->track_e > HOLD_E && e->have_last) {
        w0 = e->last_w0;
    }
    e->prev_p = 6.2831853f / w0;
    e->last_w0 = w0;
    e->have_last = 1;
    return w0;
}

/* ---- harmonics: voicing and amplitudes (py Analyzer.harmonics, IMBE's bands) ---- */

/* IMBE's w0 of b0 in 1024-point bins, Q16: w0 / 2 pi x 1024 x 65536 */
static uint32_t bins_q16(int b0)
{
    return ImbeW0q[b0] >> 6;
}

/* amplitude of harmonic l as voiced (the fitted sinusoid) or as unvoiced (the band's power), log2 Q12 */
static int32_t log2_amp(const fit_t *f, int l, int ex, int voiced)
{
    int32_t v;
    if (voiced) {
        if (f->proj[l] == 0 || f->ww[l] == 0) {
            return LOG2_FLOOR;
        }
        v = (log2_q12(f->proj[l]) - log2_q12(f->ww[l])) / 2 + (ex << 12) + LOG2_AMP_V;
    } else {
        if (f->en[l] == 0) {
            return LOG2_FLOOR;
        }
        v = log2_q12(f->en[l]) / 2 + (ex << 12) + LOG2_AMP_U;
    }
    return v > LOG2_FLOOR ? v : LOG2_FLOOR;
}

/* -> b0 of w0; voiced[1..L], log2a[1..L] (Q12) */
static int harmonics(imbe_encoder *e, const int16_t *sr, const int16_t *si, int ex, float w0, float level,
                     uint8_t *voiced, int32_t *log2a)
{
    fit_t *f = &fit;
    int b0 = imbe_b0_of((uint32_t)(w0 * 683565275.6f)); /* 2^32 / 2 pi */
    float f0 = (float)ImbeW0q[b0] * 2.3283064e-10f;   /* cycles a sample */
    int L = ImbeL[b0];
    int K = L < 37 ? (L + 2) / 3 : 12;
    harmonic_fit(sr, si, bins_q16(b0) >> FSH, L, f);
    /* the fit's error a band: IMBE's bands are three harmonics each, the last takes the rest */
    uint64_t e_b[12] = {0}, n_b[12] = {0};
    for (int l = 1; l <= L; l++) {
        int b = (l - 1) / 3 < K - 1 ? (l - 1) / 3 : K - 1;
        e_b[b] += f->en[l] - f->proj[l];
        n_b[b] += f->en[l];
    }
    float te = e->track_e * (1.0f / 65536);
    float quiet = (level - QUIET_NONE) / (QUIET_FULL - QUIET_NONE);
    quiet = quiet < 0 ? 0 : quiet > 1 ? 1 : quiet;
    if (te > (1 / f0 >= LONG_P ? LONG_MAX_E : VOICE_MAX_E)) {
        quiet = 0;
    }
    float r = (VOICE_E_NONE - te) / (VOICE_E_NONE - VOICE_E_FULL);
    quiet *= r < 0 ? 0 : r > 1 ? 1 : r;
    /* a band's decision with hysteresis against the band at the same frequency in the last frame (its centre, in
       the 500 Hz steps the AMBE encoder's state is kept in) */
    int vb[12];
    for (int b = 0; b < K; b++) {
        int c = (int)((3 * b + 2) * f0 * 16); /* 500 Hz steps */
        c = c > 15 ? 15 : c;
        float thr = VOICED_ERR * quiet * (e->prev_bands[c] ? VOICING_HYSTERESIS : 1 / VOICING_HYSTERESIS);
        vb[b] = (float)e_b[b] < thr * (float)n_b[b];
    }
    memset(e->prev_bands, 0, sizeof(e->prev_bands));
    for (int b = 0; b < K; b++) {
        int c = (int)((3 * b + 2) * f0 * 16);
        e->prev_bands[c > 15 ? 15 : c] = vb[b];
    }
    for (int l = 1; l <= L; l++) {
        int b = (l - 1) / 3 < K - 1 ? (l - 1) / 3 : K - 1;
        voiced[l] = (uint8_t)(vb[b] && f->en[l] > 0);
        log2a[l] = log2_amp(f, l, ex, voiced[l]);
    }
    return b0;
}

/* not speech: noise, all unvoiced, at a fixed fundamental (P = 32) */
static int noise_amplitudes(const int16_t *sr, const int16_t *si, int ex, uint8_t *voiced, int32_t *log2a)
{
    int b0 = imbe_b0_of(1u << 27); /* 2 pi / 32 */
    int L = ImbeL[b0];
    harmonic_fit(sr, si, bins_q16(b0) >> FSH, L, &fit);
    for (int l = 1; l <= L; l++) {
        voiced[l] = 0;
        log2a[l] = log2_amp(&fit, l, ex, 0);
    }
    return b0;
}

/* ---- the stream ---- */

void imbe_enc_init(imbe_encoder *e)
{
    memset(e, 0, sizeof(*e));
    e->E_tag[0] = e->E_tag[1] = e->E_tag[2] = -100;
    e->hist_i[0] = e->hist_i[1] = 160; /* P = 100 */
    imbe_quant_init(&e->q);
}

/* The input's DC taken out (a microphone's bias, 2500 on the DM-32: with the click at the start of a transmission it
   threw the pitch tracker off to P = 20 and every band went unvoiced; DVSI's encoder high-passes its input too):
   y = x - x1 + (1 - 1/128) y1, about 10 Hz, in integers so that all three encoders agree (py dc_block()) */
#define DC_POLE 32512 /* (1 - 1/128) in Q15 */

static int16_t dc_block(imbe_encoder *e, int16_t x)
{
    e->dc_y = ((int64_t)(x - e->dc_x) << 15) + ((e->dc_y * DC_POLE) >> 15);
    e->dc_x = x;
    int64_t y = e->dc_y >> 15;
    return (int16_t)(y > 32767 ? 32767 : y < -32768 ? -32768 : y);
}

int imbe_enc_frame(imbe_encoder *e, const int16_t in[IMBE_N], uint8_t u88[88])
{
    memmove(e->x, e->x + N, sizeof(int16_t) * (IMBE_ENC_HIST - N));
    if (e->n_in == 0) {
        e->dc_x = in[0]; /* the DC there from the start: no step into the blocker */
    }
    for (int i = 0; i < N; i++) {
        e->x[IMBE_ENC_HIST - N + i] = dc_block(e, in[i]);
    }
    e->n_in += N;
    long k = e->n_in / N - 1 - IMBE_ENC_DELAY;
    if (k < 0) {
        return 0;
    }
    /* voice activity (py vad, peak_level) */
    float edb = frame_energy_db(e, k);
    if (!e->have_level) {
        e->floor_db = e->peak_db = edb;
        e->have_level = 1;
    } else {
        e->floor_db = edb < e->floor_db + FLOOR_RISE ? edb : e->floor_db + FLOOR_RISE;
        e->peak_db = edb > e->peak_db - PEAK_DECAY ? edb : e->peak_db - PEAK_DECAY;
    }
    float thr = e->floor_db + VAD_MARGIN;
    thr = e->peak_db - VAD_BELOW_PEAK > thr ? e->peak_db - VAD_BELOW_PEAK : thr;
    thr = VAD_MIN_DB > thr ? VAD_MIN_DB : thr;
    int raw = edb > thr;
    int speech = raw || e->last_raw;
    e->last_raw = raw;
    float level = edb - e->peak_db;

    static int16_t sr[NBINS], si[NBINS];
    int ex;
    ex = spectrum(e, k, sr, si);
    int32_t log2a[58];
    uint8_t voiced[58];
    int b0;
    if (speech) {
        float w0 = pitch(e, k, sr, si);
        b0 = harmonics(e, sr, si, ex, w0, level, voiced, log2a);
    } else {
        b0 = noise_amplitudes(sr, si, ex, voiced, log2a);
    }
    imbe_quantize(&e->q, b0, voiced, log2a, u88);
    return 1;
}

int imbe_enc_frame18(imbe_encoder *e, const int16_t in[IMBE_N], uint8_t out18[18])
{
    uint8_t u[88];
    if (!imbe_enc_frame(e, in, u)) {
        return 0;
    }
    imbe_encode_frame(u, out18);
    return 1;
}
