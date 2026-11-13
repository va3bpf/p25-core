/* Tests of c/p25/p25_fec.c: random words with random errors up to and beyond what each code corrects. */
#include <stdio.h>
#include <string.h>

#include "p25_fec.h"

static uint32_t rng = 7;
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

/* a random set of `count` distinct positions below n, as a bit mask */
static uint64_t random_positions(int n, int count)
{
    uint64_t m = 0;
    while (__builtin_popcountll(m) < count) {
        m |= 1ull << (rnd() % (uint32_t)n);
    }
    return m;
}

static void test_bch(void)
{
    printf("BCH(63,16)\n");
    long beyond = 0, beyond_caught = 0;
    for (int t = 0; t < 20000; t++) {
        uint16_t d = (uint16_t)rnd(), got;
        uint64_t c = p25_bch63_encode(d);
        int ne = t % 16;
        uint64_t e = random_positions(63, ne);
        int r = p25_bch63_decode(c ^ e, &got);
        if (ne <= 11) {
            CHECK(r == ne && got == d, "%d errors: returned %d, data %04X for %04X", ne, r, got, d);
        } else {
            beyond++;
            beyond_caught += r < 0;
        }
    }
    printf("   0..11 errors corrected; 12..15 errors: %ld of %ld reported uncorrectable\n", beyond_caught, beyond);
}

static void test_golay(void)
{
    printf("Golay(24,12), Golay(18,6)\n");
    for (int t = 0; t < 20000; t++) {
        uint32_t d = rnd() & 0xFFF, got;
        int ne = t % 5;
        uint32_t c = p25_golay24_encode(d) ^ (uint32_t)random_positions(24, ne);
        int r = p25_golay24_decode(c, &got);
        if (ne <= 3) {
            CHECK(r == ne && got == d, "golay24 %d errors: %d", ne, r);
        } else {
            CHECK(r == -1, "golay24 4 errors not detected (%d)", r);
        }
        uint32_t d6 = rnd() & 0x3F;
        c = p25_golay18_encode(d6) ^ (uint32_t)random_positions(18, ne < 4 ? ne : 3);
        r = p25_golay18_decode(c, &got);
        CHECK(r == (ne < 4 ? ne : 3) && got == d6, "golay18 %d errors: %d", ne, r);
    }
    printf("   up to 3 errors corrected, 4 detected\n");
}

static void test_hamming(void)
{
    printf("Hamming(10,6)\n");
    for (uint32_t d = 0; d < 64; d++) {
        uint32_t c = p25_hamming10_encode(d), got;
        CHECK(p25_hamming10_decode(c, &got) == 0 && got == d, "clean %02X", d);
        for (int p = 0; p < 10; p++) {
            CHECK(p25_hamming10_decode(c ^ (1u << p), &got) == 1 && got == d, "%02X bit %d", d, p);
        }
    }
    printf("   every single-bit error corrected\n");
}

static void test_rs(int n, int k)
{
    int t = (n - k) / 2;
    long beyond = 0, caught = 0;
    for (int it = 0; it < 20000; it++) {
        uint8_t s[36], c[36];
        for (int i = 0; i < k; i++) {
            s[i] = (uint8_t)(rnd() & 63);
        }
        p25_rs_encode(s, n, k);
        memcpy(c, s, (size_t)n);
        int ne = it % (t + 3);
        uint64_t pos = random_positions(n, ne);
        for (int i = 0; i < n; i++) {
            if ((pos >> i) & 1) {
                c[i] ^= (uint8_t)(1 + rnd() % 63);
            }
        }
        int r = p25_rs_decode(c, n, k);
        if (ne <= t) {
            CHECK(r == ne && memcmp(c, s, (size_t)n) == 0, "RS(%d,%d) %d errors: %d", n, k, ne, r);
        } else {
            beyond++;
            caught += r < 0;
        }
    }
    printf("RS(%d,%d): up to %d symbol errors corrected; beyond: %ld of %ld reported uncorrectable\n", n, k, t, caught,
           beyond);
}

/* errors and erasures: every split with 2 e + ne <= n - k corrected, the erased symbols given as anything */
static void test_rs_erasures(int n, int k)
{
    int np = n - k;
    long tried = 0;
    for (int it = 0; it < 20000; it++) {
        uint8_t s[36], c[36], era[36];
        for (int i = 0; i < k; i++) {
            s[i] = (uint8_t)(rnd() & 63);
        }
        p25_rs_encode(s, n, k);
        memcpy(c, s, (size_t)n);
        int ne = (int)(rnd() % (np + 1));
        int nerr = (np - ne) / 2;
        nerr = nerr ? (int)(rnd() % (nerr + 1)) : 0;
        uint64_t pos = random_positions(n, ne + nerr);
        int a = 0, b = 0;
        for (int i = 0; i < n; i++) {
            if (!((pos >> i) & 1)) {
                continue;
            }
            if (a < ne) {
                era[a++] = (uint8_t)i;
                c[i] = (uint8_t)(rnd() & 63); /* maybe right by chance */
            } else {
                c[i] ^= (uint8_t)(1 + rnd() % 63);
                b++;
            }
        }
        int r = p25_rs_decode_erasures(c, n, k, era, ne);
        CHECK(r >= 0 && memcmp(c, s, (size_t)n) == 0, "RS(%d,%d) %d erasures %d errors: %d", n, k, ne, b, r);
        tried++;
    }
    printf("RS(%d,%d): %ld words with erasures and errors (2 e + ne <= %d) corrected\n", n, k, tried, np);
}

int main(void)
{
    test_bch();
    test_golay();
    test_hamming();
    test_rs(24, 12);
    test_rs(24, 16);
    test_rs(36, 20);
    test_rs_erasures(24, 12);
    test_rs_erasures(24, 16);
    test_rs_erasures(36, 20);
    printf(failures ? "%d FAILURES\n" : "all passed\n", failures);
    return failures != 0;
}
