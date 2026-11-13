/*
 * The P25 Phase 1 error-correcting codes (p25_fec.h), our own: GF(64) arithmetic, Berlekamp-Massey and Chien search
 * shared by the BCH and the Reed-Solomon decoders, Forney for the Reed-Solomon error values.
 */
#include <string.h>

#include "imbe.h"
#include "p25_fec.h"

/* ---- GF(64), x^6 + x + 1 ---- */

static uint8_t gf_exp[126], gf_log[64];
static int gf_ready;

static void gf_init(void)
{
    if (gf_ready) {
        return;
    }
    unsigned x = 1;
    for (int i = 0; i < 63; i++) {
        gf_exp[i] = (uint8_t)x;
        gf_log[x] = (uint8_t)i;
        x <<= 1;
        if (x & 64) {
            x ^= 0x43;
        }
    }
    for (int i = 63; i < 126; i++) {
        gf_exp[i] = gf_exp[i - 63];
    }
    gf_ready = 1;
}

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    return (a && b) ? gf_exp[gf_log[a] + gf_log[b]] : 0;
}

static uint8_t gf_div(uint8_t a, uint8_t b) /* b != 0 */
{
    return a ? gf_exp[gf_log[a] + 63 - gf_log[b]] : 0;
}

static uint8_t gf_pow_alpha(int e) /* alpha^e, any e */
{
    e %= 63;
    return gf_exp[e < 0 ? e + 63 : e];
}

/* Berlekamp-Massey: syndromes S[0..ns-1] (S_1..S_ns) -> error locator C (C[0] = 1); returns its degree. With ne
   erasures C comes in as their locator (degree ne), and goes out locating the errors and the erasures together. */
static int bm(const uint8_t *S, int ns, uint8_t *C, int ne)
{
    uint8_t B[24], T[24];
    memcpy(B, C, 24);
    int L = ne, m = 1;
    uint8_t b = 1;
    for (int n = ne; n < ns; n++) {
        uint8_t d = S[n];
        for (int i = 1; i <= L && i <= n; i++) {
            d ^= gf_mul(C[i], S[n - i]);
        }
        if (d == 0) {
            m++;
            continue;
        }
        uint8_t coef = gf_div(d, b);
        memcpy(T, C, 24);
        for (int i = 0; i + m < 24; i++) {
            C[i + m] ^= gf_mul(coef, B[i]);
        }
        if (2 * L <= n + ne) {
            L = n + 1 + ne - L;
            memcpy(B, T, 24);
            b = d;
            m = 1;
        } else {
            m++;
        }
    }
    return L;
}

static int berlekamp_massey(const uint8_t *S, int ns, uint8_t *C)
{
    memset(C, 0, 24);
    C[0] = 1;
    return bm(S, ns, C, 0);
}

static int chien(const uint8_t *C, int L, int *deg)
{
    int found = 0;
    for (int d = 0; d < 63; d++) {
        uint8_t v = 0;
        for (int i = 0; i <= L; i++) {
            v ^= gf_mul(C[i], gf_pow_alpha(-d * i));
        }
        if (v == 0) {
            if (found < 23) {
                deg[found] = d;
            }
            found++;
        }
    }
    return found;
}

/* ---- BCH(63,16,23) ---- */

static uint64_t bch_gen(void)
{
    /* 6331141367235453 octal: 48 bits, degree 47 */
    static const char oct[] = "6331141367235453";
    uint64_t g = 0;
    for (int i = 0; oct[i]; i++) {
        g = (g << 3) | (uint64_t)(oct[i] - '0');
    }
    return g;
}

uint64_t p25_bch63_encode(uint16_t data)
{
    uint64_t g = bch_gen(), r = (uint64_t)data << 47;
    for (int b = 62; b >= 47; b--) {
        if ((r >> b) & 1) {
            r ^= g << (b - 47);
        }
    }
    return ((uint64_t)data << 47) | r;
}

int p25_bch63_decode(uint64_t cw, uint16_t *data)
{
    gf_init();
    cw &= (1ull << 63) - 1;
    uint8_t S[22];
    int any = 0;
    for (int j = 1; j <= 22; j++) {
        uint8_t s = 0;
        for (int d = 0; d < 63; d++) {
            if ((cw >> d) & 1) {
                s ^= gf_pow_alpha(j * d);
            }
        }
        S[j - 1] = s;
        any |= s;
    }
    if (!any) {
        *data = (uint16_t)(cw >> 47);
        return 0;
    }
    uint8_t C[24];
    int L = berlekamp_massey(S, 22, C);
    int deg[23];
    if (L > 11 || chien(C, L, deg) != L) {
        *data = (uint16_t)(cw >> 47);
        return -1;
    }
    for (int i = 0; i < L; i++) {
        cw ^= 1ull << deg[i];
    }
    if (p25_bch63_encode((uint16_t)(cw >> 47)) != cw) {
        return -1;
    }
    *data = (uint16_t)(cw >> 47);
    return L;
}

/* ---- Golay(24,12), Golay(18,6) ---- */

uint32_t p25_golay24_encode(uint32_t d12)
{
    uint32_t c = imbe_golay23_encode(d12);
    return (c << 1) | (uint32_t)(__builtin_popcount(c) & 1);
}

int p25_golay24_decode(uint32_t c24, uint32_t *d12)
{
    int e;
    *d12 = imbe_golay23_decode(c24 >> 1, &e);
    int odd = __builtin_popcount(imbe_golay23_encode(*d12)) & 1, pbit = (int)(c24 & 1);
    if (odd != pbit) { /* the parity bit disagrees with the corrected word */
        if (e == 3) {
            return -1; /* 4 errors: the correction is wrong */
        }
        e++;
    }
    return e;
}

uint32_t p25_golay18_encode(uint32_t d6)
{
    return p25_golay24_encode(d6 & 0x3F) & 0x3FFFF;
}

int p25_golay18_decode(uint32_t c18, uint32_t *d6)
{
    uint32_t d12;
    int e = p25_golay24_decode(c18 & 0x3FFFF, &d12);
    if (e < 0 || (d12 & 0xFC0)) { /* the shortened bits must stay zero */
        *d6 = (c18 >> 12) & 0x3F;
        return -1;
    }
    *d6 = d12;
    return e;
}

/* ---- Hamming(10,6,3) ---- */

static const uint16_t HAM10_H[4] = {0x398, 0x354, 0x2E2, 0x1E1}; /* bit 9 first; row 0 gives syndrome bit 3 */

static uint32_t ham10_syndrome(uint32_t c10)
{
    uint32_t s = 0;
    for (int r = 0; r < 4; r++) {
        s = (s << 1) | (uint32_t)(__builtin_popcount(c10 & HAM10_H[r]) & 1);
    }
    return s;
}

uint32_t p25_hamming10_encode(uint32_t d6)
{
    uint32_t c = (d6 & 0x3F) << 4;
    return c | ham10_syndrome(c);
}

int p25_hamming10_decode(uint32_t c10, uint32_t *d6)
{
    uint32_t s = ham10_syndrome(c10);
    *d6 = (c10 >> 4) & 0x3F;
    if (!s) {
        return 0;
    }
    for (int p = 0; p < 10; p++) {
        if (ham10_syndrome(1u << p) == s) {
            *d6 = ((c10 ^ (1u << p)) >> 4) & 0x3F;
            return 1;
        }
    }
    return -1;
}

/* ---- Reed-Solomon over GF(64), shortened: sym[t] is the coefficient of x^(n-1-t) ---- */

void p25_rs_encode(uint8_t *sym, int n, int k)
{
    gf_init();
    int np = n - k;
    uint8_t g[24], par[24];
    /* g(x) = (x + alpha)(x + alpha^2)..(x + alpha^np), g[i] the coefficient of x^i */
    memset(g, 0, sizeof(g));
    g[0] = 1;
    for (int i = 1; i <= np; i++) {
        uint8_t a = gf_pow_alpha(i);
        for (int j = i; j > 0; j--) {
            g[j] = g[j - 1] ^ gf_mul(g[j], a);
        }
        g[0] = gf_mul(g[0], a);
    }
    /* the remainder of data(x) x^np by g(x): an LFSR, par[np-1] the highest degree */
    memset(par, 0, sizeof(par));
    for (int t = 0; t < k; t++) {
        uint8_t fb = (sym[t] & 63) ^ par[np - 1];
        for (int j = np - 1; j > 0; j--) {
            par[j] = par[j - 1] ^ gf_mul(fb, g[j]);
        }
        par[0] = gf_mul(fb, g[0]);
    }
    for (int j = 0; j < np; j++) {
        sym[k + j] = par[np - 1 - j];
    }
}

int p25_rs_decode_erasures(uint8_t *sym, int n, int k, const uint8_t *era, int ne)
{
    gf_init();
    int np = n - k;
    if (ne > np) {
        return -1;
    }
    uint8_t S[24];
    int any = 0;
    for (int j = 1; j <= np; j++) {
        uint8_t s = 0;
        for (int t = 0; t < n; t++) {
            s = gf_mul(s, gf_pow_alpha(j)) ^ (sym[t] & 63); /* Horner, highest degree first */
        }
        S[j - 1] = s;
        any |= s;
    }
    if (!any) {
        return 0;
    }
    uint8_t C[24]; /* the erasures' locator: the product of (1 + alpha^(n-1-era[i]) x) */
    memset(C, 0, sizeof(C));
    C[0] = 1;
    for (int i = 0; i < ne; i++) {
        uint8_t x = gf_pow_alpha(n - 1 - era[i]);
        for (int j = i + 1; j > 0; j--) {
            C[j] ^= gf_mul(x, C[j - 1]);
        }
    }
    int L = bm(S, np, C, ne);
    int deg[23];
    if (2 * L > np + ne || chien(C, L, deg) != L) {
        return -1;
    }
    for (int i = 0; i < L; i++) {
        if (deg[i] >= n) {
            return -1; /* an error in the shortened (zero) part: too many errors */
        }
    }
    /* Forney: Omega(x) = S(x) C(x) mod x^np, S(x) = sum S_j x^(j-1); e = Omega(X^-1) / C'(X^-1) */
    uint8_t om[24];
    memset(om, 0, sizeof(om));
    for (int i = 0; i < np; i++) {
        for (int j = 0; j <= L && i + j < np; j++) {
            om[i + j] ^= gf_mul(S[i], C[j]);
        }
    }
    uint8_t fix[36];
    memset(fix, 0, sizeof(fix));
    for (int i = 0; i < L; i++) {
        uint8_t xinv = gf_pow_alpha(-deg[i]), num = 0, den = 0, p = 1;
        for (int j = 0; j < np; j++) {
            num ^= gf_mul(om[j], p);
            p = gf_mul(p, xinv);
        }
        p = 1;
        for (int j = 1; j <= L; j += 2) { /* C'(x): the odd terms, x^(j-1) */
            den ^= gf_mul(C[j], p);
            p = gf_mul(p, gf_mul(xinv, xinv));
        }
        if (den == 0) {
            return -1;
        }
        fix[n - 1 - deg[i]] = gf_div(num, den);
    }
    for (int t = 0; t < n; t++) {
        sym[t] ^= fix[t];
    }
    return L;
}

int p25_rs_decode(uint8_t *sym, int n, int k)
{
    return p25_rs_decode_erasures(sym, n, k, NULL, 0);
}
