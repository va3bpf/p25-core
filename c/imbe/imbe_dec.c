/*
 * IMBE decoder: 18 bytes -> 160 samples. The synthesis is dmr-vocoder's (c/ambe/ambe_dec.c, our own): full-strength
 * spectral enhancement, voiced harmonics as sinusoids (phase-continuous below the 8th harmonic, a fixed phase offset
 * each), unvoiced bands as shaped noise, a 4th-order high-pass; the unvoiced gain and the high-pass retuned for IMBE
 * against imbe_vocoder (below). Frame repeat and muting as mbelib: a frame with
 * an invalid pitch or more than 5 corrected bits repeats the last one, the fourth repeat in a row mutes.
 *
 * Integer: phases in uint32 (2^32 = 2 pi), the cosine from a Q12 table, the output in Q4 until the high-pass (Q8);
 * the amplitudes Ml in Q6, the enhancement and the noise gains in log2 (Q12).
 */
#include <string.h>

#include "imbe.h"
#include "imbe_fx.h"

#define N IMBE_N
#define INTERP_MAX 8
#define NM 512               /* the noise's IFFT */
#define AMAX ((1 << 19) - 1) /* largest amplitude, Q4 (32768 PCM): x Q12 cosine stays in 32 bits */
/* the unvoiced noise's gain, log2 in Q12: 1.7 (dmr-vocoder's 0.77 for AMBE left the top band 2.6 dB below the reference
   decoder's, its sound dull; py/speech_test.py PESQ at a 120 Hz high-pass: 0.77 2.824, 1.0 2.870, 1.3 2.899, 1.5 2.911,
   1.7 2.917, 2.0 2.907, 2.2 2.902) */
#ifndef IMBE_UV_LOG2_Q12
#define IMBE_UV_LOG2_Q12 3136
#endif
#define REPEAT_ERRORS 5      /* more corrected bits than this: repeat the last frame (mbelib) */
#define MUTE_REPEATS 3       /* more repeats in a row than this: mute */

/* the high-pass: 4th-order Butterworth as two RBJ sections (Q 0.5412, 1.3066): b0 b1 b2 a1 a2, Q14 */
/* 250 Hz (dmr-vocoder's 120 left IMBE's output 4 dB heavier below 250 Hz than the reference decoder's; py/speech_test.py
   PESQ with UV gain 1.7: 120 Hz 2.917, 200 Hz 2.940, 250 Hz 2.960, 300 Hz 2.981 but 10 dB less below 250 Hz) */
#ifndef IMBE_HP_HZ
#define IMBE_HP_HZ 250
#endif
#if IMBE_HP_HZ == 120
static const int16_t HP[2][5] = {{15040, -30080, 15040, -30013, 13763}, {15779, -31559, 15779, -31489, 15245}};
#elif IMBE_HP_HZ == 150
static const int16_t HP[2][5] = {{14728, -29456, 14728, -29353, 13174}, {15624, -31249, 15624, -31140, 14974}};
#elif IMBE_HP_HZ == 200
static const int16_t HP[2][5] = {{14227, -28454, 14227, -28278, 12246}, {15363, -30727, 15363, -30537, 14533}};
#elif IMBE_HP_HZ == 250
static const int16_t HP[2][5] = {{13749, -27497, 13749, -27230, 11380}, {15099, -30199, 15099, -29906, 14108}};
#elif IMBE_HP_HZ == 300
static const int16_t HP[2][5] = {{13291, -26582, 13291, -26210, 10571}, {14833, -29665, 14833, -29250, 13697}};
#else
#error "IMBE_HP_HZ: 120, 150, 200, 250 or 300"
#endif

int32_t imbe_scratch[IMBE_SCRATCH_WORDS];

static uint32_t rnd(imbe_decoder *d)
{
    uint32_t x = d->rng; /* xorshift32 */
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    d->rng = x;
    return x;
}

/* about normal, sd 1180: the sum of four uniform numbers */
static int32_t rnd_normal(imbe_decoder *d)
{
    int32_t s = 0;
    for (int i = 0; i < 4; i++) {
        s += (int32_t)(rnd(d) >> 16) - 32768;
    }
    return s >> 5;
}

/* (a * b) >> 12 for a 32-bit a and a Q12 b */
static inline int32_t mulq12(int32_t a, int32_t b)
{
    return (a >> 12) * b + (((a & 0xFFF) * b) >> 12);
}

static inline int32_t mulq14(int32_t a, int32_t b)
{
    return (a >> 14) * b + (((a & 0x3FFF) * b) >> 14);
}

/* Ml (Q6) -> the output amplitude GAIN x Ml, Q4 */
static int32_t amp_q4(int32_t ml)
{
    int32_t a = mulq12(ml, 2406); /* 2.35 / 4 */
    return a > AMAX ? AMAX : a;
}

/* ---- spectral amplitude enhancement (mbe_spectralAmpEnhance, full strength), integer ----
   W_l^4 = M_l^2 0.96 pi (1 + rho^2 - 2 rho cos(w0 l)) / (w0 Rm0 (1 - rho^2)), rho = Rm1 / Rm0: mbelib's formula
   divided through by Rm0^3, evaluated in log2 (Q12) */

static void spectral_amp_enhance(imbe_parms *cur)
{
    int L = cur->L;
    int32_t mx = 0;
    for (int l = 1; l <= L; l++) {
        mx = cur->Ml[l] > mx ? cur->Ml[l] : mx;
    }
    if (mx == 0) {
        return;
    }
    int sh = (31 - __builtin_clz((uint32_t)mx)) - 11; /* m_l = Ml >> sh < 2^12: the sums fit 32 bits */
    int32_t m[58];
    uint32_t Rm0 = 0;
    int32_t Rm1 = 0;
    for (int l = 1; l <= L; l++) {
        m[l] = sh >= 0 ? cur->Ml[l] >> sh : cur->Ml[l] << -sh;
        int32_t m2 = m[l] * m[l];
        Rm0 += (uint32_t)m2;
        Rm1 += (m2 >> 12) * cos_q12(cur->w0q * (uint32_t)l);
    }
    if (Rm0 == 0) {
        return;
    }
    int32_t rho = (int32_t)(((int64_t)Rm1 << 15) / (int64_t)Rm0); /* Q15 */
    int32_t rho2 = (rho * rho) >> 15, one_m = 32768 - rho2;
    if (one_m <= 0) {
        return;
    }
    /* log2 of the factors common to every l: 0.96 pi / w0 = 0.48 x 2^32 / w0q; / Rm0; / (1 - rho^2) */
    int32_t common = 126735 - log2_q12(cur->w0q) - log2_q12(Rm0) - (log2_q12((uint32_t)one_m) - (15 << 12));
    uint32_t s = 0;
    int32_t f[58];
    for (int l = 1; l <= L; l++) {
        f[l] = 1 << 14;
        if (m[l] != 0 && 8 * l > L) {
            int32_t c = cos_q12(cur->w0q * (uint32_t)l);
            int32_t a = 32768 + rho2 - ((2 * rho * c) >> 12); /* 1 + rho^2 - 2 rho cos, Q15 */
            int32_t lw = -4096;                              /* log2 0.5 */
            if (a > 0) {
                lw = (log2_q12((uint32_t)(m[l] * m[l])) + log2_q12((uint32_t)a) - (15 << 12) + common) >> 2;
                lw = lw > 1077 ? 1077 : lw < -4096 ? -4096 : lw; /* 0.5 .. 1.2 */
            }
            f[l] = (int32_t)exp2_q12(lw, 14);
        }
        int32_t nm = (m[l] * f[l]) >> 14;
        s += (uint32_t)(nm * nm);
    }
    /* energy kept: g = sqrt(Rm0 / s) */
    int32_t g = s ? (int32_t)exp2_q12((log2_q12(Rm0) - log2_q12(s)) / 2, 14) : (1 << 14);
    for (int l = 1; l <= L; l++) {
        int32_t k = (f[l] * g) >> 14;
        cur->Ml[l] = mulq14(cur->Ml[l], k > 32767 ? 32767 : k);
    }
}

/* ---- voiced synthesis, into out (Q4) ---- */

static void synth_voiced(imbe_decoder *d, imbe_parms *cur, imbe_parms *prev, int32_t *out)
{
    int num_uv = 0;
    for (int l = 1; l <= cur->L; l++) {
        num_uv += cur->Vl[l] == 0;
    }
    uint32_t cw0 = cur->w0q, pw0 = prev->w0q;
    int maxl;
    if (cur->L > prev->L) {
        maxl = cur->L;
        for (int l = prev->L + 1; l <= maxl; l++) {
            prev->Ml[l] = 0;
            prev->Vl[l] = 1;
        }
    } else {
        maxl = prev->L;
        for (int l = cur->L + 1; l <= maxl; l++) {
            cur->Ml[l] = 0;
            cur->Vl[l] = 1;
        }
    }
    for (int l = 1; l <= 56; l++) {
        cur->PSIl[l] = prev->PSIl[l] + (pw0 + cw0) * (uint32_t)(l * N / 2);
        cur->PHIl[l] = cur->PSIl[l];
        if (l > cur->L / 4) { /* + num_uv x a random phase / L */
            cur->PHIl[l] += (uint32_t)((((int32_t)rnd(d) >> 6) * num_uv / cur->L) << 6);
        }
        cur->PHIl[l] += (uint32_t)ImbePhaseDispQ32[l];
    }
    int32_t *acc_p = imbe_scratch, *acc_c = imbe_scratch + N;
    memset(acc_p, 0, sizeof(int32_t) * 2 * N);
    int dw0_small = (int32_t)(cw0 - pw0) < (int32_t)(cw0 / 10) && (int32_t)(pw0 - cw0) < (int32_t)(cw0 / 10);
    for (int l = 1; l <= maxl; l++) {
        int32_t ap = amp_q4(prev->Ml[l]), ac = amp_q4(cur->Ml[l]);
        if (prev->Vl[l] && cur->Vl[l] && l < INTERP_MAX && dw0_small) {
            /* one sinusoid: amplitude ramp, the phase lands on PHIl at n = N (IMBE eq. 137, 138) */
            int32_t dphi = (int32_t)(cur->PHIl[l] - prev->PHIl[l] - (pw0 + cw0) * (uint32_t)(l * N / 2));
            int32_t c = (int32_t)(cw0 - pw0) * l / (2 * N); /* the phase's second difference / 2 */
            uint32_t th = prev->PHIl[l], inc = pw0 * (uint32_t)l + (uint32_t)(dphi / N) + (uint32_t)c;
            int32_t amp = ap << 8, damp = ((ac - ap) << 8) / N;
            for (int n = 0; n < N; n++) {
                out[n] += ((amp >> 8) * cos_q12(th)) >> 12;
                th += inc;
                inc += (uint32_t)(2 * c);
                amp += damp;
            }
            continue;
        }
        if (prev->Vl[l] && ap) {
            uint32_t th = prev->PHIl[l], step = pw0 * (uint32_t)l;
            for (int n = 0; n < N; n++) {
                acc_p[n] += (ap * cos_q12(th)) >> 12;
                th += step;
            }
        }
        if (cur->Vl[l] && ac) {
            uint32_t step = cw0 * (uint32_t)l, th = cur->PHIl[l] - step * N;
            for (int n = 0; n < N; n++) {
                acc_c[n] += (ac * cos_q12(th)) >> 12;
                th += step;
            }
        }
    }
    for (int n = 0; n < N; n++) {
        out[n] += mul32x16(acc_p[n], ImbeWsQ15[N + n]) + mul32x16(acc_c[n], ImbeWsQ15[n]);
    }
}

/* ---- unvoiced synthesis: white noise kept in the unvoiced bands, each with the power of a sinusoid of amplitude
   Ml, 2N samples a frame (of an NM-point IFFT), cross-faded with sqrt(Ws) ---- */

/* in-place inverse FFT (e^+i), no scaling, Q12 twiddles from the cosine table */
static void ifft(int32_t *re, int32_t *im, int n)
{
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
        uint32_t step = 0xFFFFFFFFu / (uint32_t)len + 1; /* 2^32 / len */
        for (int j = 0; j < len / 2; j++) {
            uint32_t ph = step * (uint32_t)j;
            int32_t wr = cos_q12(ph), wi = cos_q12(ph - (1u << 30)); /* e^(+i 2 pi j / len) */
            for (int i = j; i < n; i += len) {
                int32_t *ar = &re[i], *ai = &im[i], *br = &re[i + len / 2], *bi = &im[i + len / 2];
                int32_t vr = mulq12(*br, wr) - mulq12(*bi, wi);
                int32_t vi = mulq12(*br, wi) + mulq12(*bi, wr);
                *br = *ar - vr;
                *bi = *ai - vi;
                *ar += vr;
                *ai += vi;
            }
        }
    }
}

static void synth_noise(imbe_decoder *d, const imbe_parms *p, int32_t *out)
{
    int32_t *xr = imbe_scratch, *xi = imbe_scratch + NM; /* after synth_voiced */
    uint8_t in[NM / 2];
    memset(xr, 0, sizeof(int32_t) * 2 * NM);
    memset(in, 0, sizeof(in));
    int any = 0;
    uint32_t wb = p->w0q >> 7; /* w0 in NM-point bins, Q16 */
    for (int k = 1; k < NM / 2; k++) {
        xr[k] = rnd_normal(d);
        xi[k] = rnd_normal(d);
    }
    for (int l = 1; l <= p->L; l++) {
        if (p->Vl[l] || p->Ml[l] == 0) {
            continue;
        }
        /* bins with (l - 1/2) w0 <= w_k < (l + 1/2) w0, 0 < k < NM / 2 */
        uint32_t lo = (uint32_t)(2 * l - 1) * wb / 2, hi = (uint32_t)(2 * l + 1) * wb / 2; /* Q16 */
        int k0 = (int)((lo + 0xFFFF) >> 16);                                               /* k >= lo */
        if (k0 < 1) {
            k0 = 1;
        }
        int k1 = (int)((hi - 1) >> 16); /* k < hi */
        if (k1 > NM / 2 - 1) {
            k1 = NM / 2 - 1;
        }
        uint32_t e = 0;
        for (int k = k0; k <= k1; k++) {
            e += (uint32_t)(xr[k] * xr[k] + xi[k] * xi[k]);
        }
        if (e == 0) {
            continue;
        }
        /* u = 2 / NM sum Re(S_k e^(...)), S_k = X_k Ml NM / (2 sqrt(e)): the scale per X_k is
           UV_GAIN GAIN Ml / sqrt(e) in Q4, Ml in Q6: log2(UV_GAIN) + log2(2.35 / 4) (-3144) */
        int32_t lg = log2_q12((uint32_t)p->Ml[l]) + IMBE_UV_LOG2_Q12 - 3144 - log2_q12(e) / 2;
        uint32_t g = exp2_q12(lg, 8);
        g = g > (1u << 18) ? (1u << 18) : g;
        for (int k = k0; k <= k1; k++) {
            xr[k] = (xr[k] * (int32_t)g) >> 8;
            xi[k] = (xi[k] * (int32_t)g) >> 8;
            in[k] = 1;
        }
        any = 1;
    }
    for (int k = 1; k < NM / 2; k++) { /* bins outside the unvoiced bands carry nothing */
        if (!in[k]) {
            xr[k] = xi[k] = 0;
        }
    }
    if (any) {
        ifft(xr, xi, NM);
    } else {
        memset(xr, 0, sizeof(int32_t) * NM);
    }
    for (int n = 0; n < N; n++) {
        out[n] += mul32x16(d->noise_tail[n], ImbeSqrtWsQ15[N + n]) + mul32x16(xr[n], ImbeSqrtWsQ15[n]);
        d->noise_tail[n] = xr[N + n];
    }
}

/* ---- high-pass, Q8 ---- */

static void highpass(imbe_decoder *d, int32_t *x)
{
    for (int s = 0; s < 2; s++) {
        const int16_t *c = HP[s];
        int32_t *st = d->hp[s];
        for (int n = 0; n < N; n++) {
            int32_t v = x[n];
            int32_t o = mulq14(v, c[0]) + mulq14(st[0], c[1]) + mulq14(st[1], c[2]) - mulq14(st[2], c[3]) -
                        mulq14(st[3], c[4]);
            st[1] = st[0];
            st[0] = v;
            st[3] = st[2];
            st[2] = o;
            x[n] = o;
        }
    }
}

static void reset_voice(imbe_decoder *d)
{
    imbe_parms_init(&d->cur);
    imbe_parms_init(&d->prev);
    imbe_parms_init(&d->prev_enh);
    memset(d->noise_tail, 0, sizeof(d->noise_tail));
}

void imbe_dec_init(imbe_decoder *d)
{
    memset(d, 0, sizeof(*d));
    d->rng = 0x12345678;
    reset_voice(d);
}

int imbe_dec_parms_frame(imbe_decoder *d, const uint8_t u88[88], const imbe_errors *err, int16_t pcm[IMBE_N])
{
    int32_t out[N];
    int kind;
    memset(out, 0, sizeof(out));
    imbe_parms *cur = &d->cur, *prev = &d->prev;
    int bad = imbe_decode_parms(u88, cur, prev);
    if (bad || err->etotal > REPEAT_ERRORS) {
        *cur = *prev;
        cur->repeat++;
    } else {
        cur->repeat = 0;
    }
    if (cur->repeat <= MUTE_REPEATS) {
        *prev = *cur;
        spectral_amp_enhance(cur);
        synth_voiced(d, cur, &d->prev_enh, out);
        synth_noise(d, cur, out);
        d->prev_enh = *cur;
        kind = cur->repeat ? 'R' : 'V';
    } else {
        reset_voice(d);
        kind = 'M';
    }
    for (int n = 0; n < N; n++) {
        out[n] <<= 4; /* Q4 -> Q8 */
    }
    highpass(d, out);
    for (int n = 0; n < N; n++) {
        int32_t v = (out[n] + 128) >> 8;
        pcm[n] = (int16_t)(v > 32760 ? 32760 : v < -32760 ? -32760 : v);
    }
    return kind;
}

int imbe_dec_frame(imbe_decoder *d, const uint8_t in18[18], int16_t pcm[IMBE_N])
{
    uint8_t u88[88];
    imbe_errors err;
    imbe_decode_frame(in18, u88, &err);
    return imbe_dec_parms_frame(d, u88, &err, pcm);
}
