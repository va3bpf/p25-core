/*
 * Tests of c/imbe on the PC (c/build.sh), the reference being mbelib (ref/mbelib, linked here only):
 *   1. Golay(23,12): every error pattern of up to 3 bits on random words; our parity = mbelib's golayGenerator
 *   2. Hamming(15,11): every single-bit error
 *   3. the frame layer: random parameters -> 144 bits -> random correctable errors -> the same parameters
 *   4. the frame layer against mbelib's (bit order, scrambling): the same 88 bits from clean frames
 *   5. the parameters against mbe_decodeImbe4400Parms on random frames: L, voicing exact, w0 and log2 Ml close
 *   6. the output level against mbelib's synthesis (RMS ratio)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "imbe.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Tests 4..6 link mbelib (-DWITH_MBELIB, c/build.sh mbelib); 1..3 need only our code */
#ifdef WITH_MBELIB
#include "mbelib.h"
#endif

/* mbelib's golayGenerator (ref/mbelib/ecc_const.h): the parity bits of each data bit, MSB first */
static const int golayGenerator[12] = {0x63a, 0x31d, 0x7b4, 0x3da, 0x1ed, 0x6cc,
                                       0x366, 0x1b3, 0x6e3, 0x54b, 0x49f, 0x475};

static uint32_t rng = 1;
static uint32_t rnd(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    return rng;
}

static int failures;
#define CHECK(cond, ...)                                                                                               \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            if (failures++ < 20) {                                                                                     \
                printf("  FAIL: " __VA_ARGS__);                                                                        \
                printf("\n");                                                                                          \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

static void test_golay(void)
{
    printf("1. Golay(23,12)\n");
    for (int i = 0; i < 12; i++) {
        uint32_t p = imbe_golay23_encode(1u << (11 - i)) & 0x7FF; /* data bit i, MSB first, as mbelib's table */
        CHECK(p == (uint32_t)golayGenerator[i], "parity of data bit %d: %03X, mbelib %03X", i, p,
              golayGenerator[i]);
    }
    long n = 0;
    for (int t = 0; t < 200; t++) {
        uint32_t d = rnd() & 0xFFF, c = imbe_golay23_encode(d);
        for (int a = -1; a < 23; a++) {
            for (int b = a; b < 23; b++) {
                for (int e = b; e < 23; e++) {
                    if ((a == b && a != -1) || (b == e && b != -1)) {
                        continue;
                    }
                    uint32_t err = (a >= 0 ? 1u << a : 0) | (b >= 0 ? 1u << b : 0) | (e >= 0 ? 1u << e : 0);
                    int k;
                    uint32_t got = imbe_golay23_decode(c ^ err, &k);
                    CHECK(got == d && k == __builtin_popcount(err), "data %03X errors %06X: %03X (%d)", d, err, got,
                          k);
                    n++;
                }
            }
        }
    }
    printf("   %ld words decoded\n", n);
}

static void test_hamming(void)
{
    printf("2. Hamming(15,11)\n");
    for (uint32_t d = 0; d < 2048; d++) {
        uint32_t c = imbe_hamming15_encode(d);
        int k;
        CHECK(imbe_hamming15_decode(c, &k) == d && k == 0, "clean %03X", d);
        for (int p = 0; p < 15; p++) {
            CHECK(imbe_hamming15_decode(c ^ (1u << p), &k) == d && k == 1, "data %03X bit %d", d, p);
        }
    }
    printf("   2048 words x 16 patterns\n");
}

static const int CW_LEN[8] = {23, 23, 23, 23, 15, 15, 15, 7};

static void random_u(uint8_t u[88])
{
    for (;;) {
        for (int i = 0; i < 88; i++) {
            u[i] = (uint8_t)(rnd() & 1);
        }
        int b0 = (u[0] << 7) | (u[1] << 6) | (u[2] << 5) | (u[3] << 4) | (u[4] << 3) | (u[5] << 2) | (u[85] << 1) |
                 u[86];
        if (b0 <= 207) {
            return;
        }
    }
}

static void flip(uint8_t *buf, int pos)
{
    buf[pos >> 3] ^= (uint8_t)(0x80 >> (pos & 7));
}

static void test_frame(void)
{
    printf("3. frame layer, correctable errors\n");
    for (int t = 0; t < 20000; t++) {
        uint8_t u[88], f[18], v[88];
        random_u(u);
        imbe_encode_frame(u, f);
        /* up to 3 errors in each of c0..c3, 1 in each of c4..c6, none in c7 */
        int pos = 0, nerr = 0;
        for (int w = 0; w < 7; w++) {
            int m = w < 4 ? (int)(rnd() % 4) : (int)(rnd() % 2);
            int used = 0;
            for (int k = 0; k < m; k++) {
                int b = (int)(rnd() % (uint32_t)CW_LEN[w]);
                if (used & (1 << b)) {
                    continue;
                }
                used |= 1 << b;
                flip(f, pos + b);
                nerr++;
            }
            pos += CW_LEN[w];
        }
        imbe_errors e;
        imbe_decode_frame(f, v, &e);
        CHECK(memcmp(u, v, 88) == 0, "frame %d not recovered", t);
        CHECK(e.etotal == nerr, "frame %d: %d errors counted, %d made", t, e.etotal, nerr);
    }
    printf("   20000 frames\n");
}

#ifdef WITH_MBELIB
/* our 18 bytes -> mbelib's imbe_fr[8][23] (imbe_fr[i][j] = bit j of c_i) */
static void to_mbelib(const uint8_t f[18], char fr[8][23])
{
    int pos = 0;
    memset(fr, 0, 8 * 23);
    for (int i = 0; i < 8; i++) {
        for (int j = CW_LEN[i] - 1; j >= 0; j--, pos++) {
            fr[i][j] = (char)((f[pos >> 3] >> (7 - (pos & 7))) & 1);
        }
    }
}

static void test_frame_mbelib(void)
{
    printf("4. frame layer against mbelib (clean frames)\n");
    for (int t = 0; t < 20000; t++) {
        uint8_t u[88], f[18], v[88];
        char fr[8][23], d[88];
        random_u(u);
        imbe_encode_frame(u, f);
        imbe_errors e;
        imbe_decode_frame(f, v, &e);
        to_mbelib(f, fr);
        int e0 = mbe_eccImbe7200x4400C0(fr);
        mbe_demodulateImbe7200x4400Data(fr);
        int e1 = mbe_eccImbe7200x4400Data(fr, d);
        int same = 1;
        for (int i = 0; i < 88; i++) {
            same &= d[i] == (char)v[i];
        }
        CHECK(same && e0 == 0 && e1 == 0, "frame %d differs from mbelib (errors %d %d)", t, e0, e1);
    }
    printf("   20000 frames\n");
}

static void test_parms(void)
{
    printf("5. parameters against mbelib\n");
    imbe_parms cur, prev;
    mbe_parms mc, mp, me;
    imbe_parms_init(&cur);
    imbe_parms_init(&prev);
    mbe_initMbeParms(&mc, &mp, &me);
    double worst = 0, sum = 0, worst_w0 = 0;
    long cnt = 0, big = 0, l_diff = 0, v_diff = 0;
    for (int t = 0; t < 20000; t++) {
        uint8_t u[88];
        char d[88];
        random_u(u);
        for (int i = 0; i < 88; i++) {
            d[i] = (char)u[i];
        }
        mp.log2Ml[0] = mp.log2Ml[1]; /* as the standard (and ours) */
        int r1 = imbe_decode_parms(u, &cur, &prev);
        int r2 = mbe_decodeImbe4400Parms(d, &mc, &mp);
        CHECK(r1 == r2, "frame %d: valid %d, mbelib %d", t, r1, r2);
        if (r1 || r2) {
            continue;
        }
        if (cur.L != mc.L) {
            l_diff++;
            continue;
        }
        double w0 = cur.w0q * (2 * M_PI / 4294967296.0);
        double dw = fabs(w0 - mc.w0) / mc.w0;
        worst_w0 = dw > worst_w0 ? dw : worst_w0;
        int frame_big = 0;
        for (int l = 1; l <= cur.L; l++) {
            v_diff += cur.Vl[l] != mc.Vl[l];
            double a = cur.log2Ml[l] / 4096.0, b = mc.log2Ml[l];
            double e = fabs(a - b);
            worst = e > worst ? e : worst;
            sum += e;
            cnt++;
            frame_big |= e > 0.05;
        }
        big += frame_big;
        prev = cur;
        mbe_moveMbeParms(&mc, &mp);
        /* restart both from the same state now and then, so one frame's difference does not run on */
        if ((t % 50) == 49) {
            for (int l = 0; l < 57; l++) {
                prev.log2Ml[l] = (int32_t)lrint(mp.log2Ml[l] * 4096);
            }
        }
    }
    printf("   L differs: %ld frames, voicing differs: %ld harmonics\n", l_diff, v_diff);
    printf("   w0: worst relative error %.2e\n", worst_w0);
    printf("   log2 Ml: mean |error| %.5f, worst %.4f, frames with an error above 0.05: %ld\n", sum / cnt, worst, big);
    CHECK(l_diff == 0 && v_diff == 0, "L / voicing differ");
    CHECK(worst_w0 < 1e-4, "w0");
}

static void test_level(void)
{
    printf("6. output level against mbelib\n");
    imbe_decoder dec;
    imbe_dec_init(&dec);
    mbe_parms mc, mp, me;
    mbe_initMbeParms(&mc, &mp, &me);
    double e_ours = 0, e_mbe = 0;
    /* a slowly varying frame sequence: the same random frame 10 times, then the next */
    uint8_t u[88];
    for (int t = 0; t < 2000; t++) {
        if (t % 10 == 0) {
            random_u(u);
        }
        imbe_errors err = {0, 0, 0};
        int16_t pcm[IMBE_N];
        short mpcm[160];
        char d[88], es[64];
        int e1 = 0, e2 = 0;
        for (int i = 0; i < 88; i++) {
            d[i] = (char)u[i];
        }
        imbe_dec_parms_frame(&dec, u, &err, pcm);
        mbe_processImbe4400Data(mpcm, &e1, &e2, es, d, &mc, &mp, &me, 3);
        for (int n = 0; n < IMBE_N; n++) {
            e_ours += (double)pcm[n] * pcm[n];
            e_mbe += (double)mpcm[n] * mpcm[n];
        }
    }
    printf("   RMS ours / mbelib: %.3f (%.1f dB)\n", sqrt(e_ours / e_mbe), 10 * log10(e_ours / e_mbe));
}
#endif

int main(void)
{
    test_golay();
    test_hamming();
    test_frame();
#ifdef WITH_MBELIB
    test_frame_mbelib();
    test_parms();
    test_level();
#else
    printf("(4..6 against mbelib: build with c/build.sh mbelib)\n");
#endif
    printf(failures ? "%d FAILURES\n" : "all passed\n", failures);
    return failures != 0;
}
