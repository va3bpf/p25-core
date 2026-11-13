/*
 * Tests of c/p25/p25_frame.c: a whole call (HDU, LDU1 / LDU2 pairs, TDULC) built by the transmit side, behind some
 * random dibits, with random bit errors, through the receiver: the data units, NAC, link control, encryption sync
 * and every IMBE frame's parameters (after c/imbe's error correction).
 */
#include <stdio.h>
#include <string.h>

#include "imbe.h"
#include "p25_frame.h"

static uint32_t rng = 11;
static uint32_t rnd(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    return rng;
}

#define PAIRS 10
#define STREAM (300 + 396 + PAIRS * 2 * 864 + 216 + 100)
#define NAC 0x293

static uint8_t stream[STREAM];
static uint8_t u_sent[PAIRS * 2 * 9][88];

static int build(void)
{
    int n = 0;
    for (int i = 0; i < 300; i++) {
        stream[n++] = (uint8_t)(rnd() & 3); /* noise before the call */
    }
    p25_es es = {{0}, P25_ALGID_CLEAR, 0, 0};
    p25_lc lc;
    memset(&lc, 0, sizeof(lc));
    lc.lco = 0x00;
    lc.tgid = 1;
    lc.src = 1234567;
    p25_lc_pack(&lc);
    n += p25_tx_hdu(stream + n, NAC, &es, 0, 1, 2);
    for (int p = 0; p < PAIRS; p++) {
        for (int half = 0; half < 2; half++) {
            uint8_t imbe[9][18];
            for (int k = 0; k < 9; k++) {
                uint8_t *u = u_sent[(p * 2 + half) * 9 + k];
                for (int i = 0; i < 88; i++) {
                    u[i] = (uint8_t)(rnd() & 1);
                }
                imbe_encode_frame(u, imbe[k]);
            }
            uint8_t lsd[2] = {(uint8_t)p, (uint8_t)(p + 1)};
            n += p25_tx_ldu(stream + n, NAC, half ? P25_DUID_LDU2 : P25_DUID_LDU1, (const uint8_t(*)[18])imbe, &lc,
                            &es, lsd, 2);
        }
    }
    n += p25_tx_tdulc(stream + n, NAC, &lc, 2);
    for (int i = 0; i < 100; i++) {
        stream[n++] = (uint8_t)(rnd() & 3);
    }
    return n;
}

typedef struct {
    int units, hdu, ldu, tdulc, lc_ok, es_ok, imbe_ok, imbe_total;
} result;

static result run(int n, double ber)
{
    static uint8_t rx_stream[STREAM];
    memcpy(rx_stream, stream, (size_t)n);
    uint32_t threshold = (uint32_t)(ber * 4294967295.0);
    for (int i = 0; i < n; i++) {
        for (int b = 0; b < 2; b++) {
            if (rnd() < threshold) {
                rx_stream[i] ^= (uint8_t)(1 << b);
            }
        }
    }
    static p25_rx rx;
    p25_rx_init(&rx, P25_NAC_ANY);
    result r;
    memset(&r, 0, sizeof(r));
    int ldu_index = 0;
    for (int i = 0; i < n; i++) {
        int duid = p25_rx_dibit(&rx, rx_stream[i]);
        if (duid < 0) {
            continue;
        }
        const p25_frame *f = &rx.frame;
        r.units++;
        if (f->nac != NAC) {
            continue;
        }
        if (duid == P25_DUID_HDU) {
            r.hdu++;
            r.es_ok += f->es.ok && f->es.algid == P25_ALGID_CLEAR && f->tgid == 1;
        } else if (duid == P25_DUID_LDU1 || duid == P25_DUID_LDU2) {
            /* which LDU of the call this is: from the LSD byte where it is intact, else the running count */
            r.ldu++;
            if (duid == P25_DUID_LDU1) {
                r.lc_ok += f->lc.ok && f->lc.tgid == 1 && f->lc.src == 1234567;
            } else {
                r.es_ok += f->es.ok && f->es.algid == P25_ALGID_CLEAR;
            }
            for (int k = 0; k < 9; k++) {
                uint8_t u[88];
                imbe_errors e;
                imbe_decode_frame(f->imbe[k], u, &e);
                r.imbe_total++;
                if (ldu_index < PAIRS * 2) {
                    r.imbe_ok += memcmp(u, u_sent[ldu_index * 9 + k], 88) == 0;
                }
            }
            ldu_index++;
        } else if (duid == P25_DUID_TDULC) {
            r.tdulc++;
            r.lc_ok += f->lc.ok && f->lc.tgid == 1 && f->lc.src == 1234567;
        }
    }
    return r;
}

/* p25_tx_ldu_ready(): the LDU's first air dibits don't change with the voice frames from there on (sent as soon as
   their frames are encoded), and not much more than that is held back */
static int check_ldu_ready(void)
{
    static uint8_t a[9][18], b[9][18], la[864], lb[864];
    p25_lc lc;
    p25_es es;
    memset(&lc, 0, sizeof(lc));
    lc.tgid = 1;
    lc.src = 3027527;
    p25_lc_pack(&lc);
    memset(&es, 0, sizeof(es));
    es.algid = P25_ALGID_CLEAR;
    int bad = 0;
    for (int duid = P25_DUID_LDU1; duid <= P25_DUID_LDU2; duid += P25_DUID_LDU2 - P25_DUID_LDU1) {
        for (int n = 0; n <= 9; n++) {
            for (int k = 0; k < 9; k++) {
                for (int i = 0; i < 18; i++) {
                    a[k][i] = (uint8_t)rnd();
                    b[k][i] = (k < n) ? a[k][i] : (uint8_t)~a[k][i];
                }
            }
            p25_tx_ldu(la, NAC, duid, (const uint8_t(*)[18])a, &lc, &es, NULL, 2);
            p25_tx_ldu(lb, NAC, duid, (const uint8_t(*)[18])b, &lc, &es, NULL, 2);
            int first = 864;
            for (int i = 0; i < 864; i++) {
                if (la[i] != lb[i]) {
                    first = i;
                    break;
                }
            }
            int ready = p25_tx_ldu_ready(n);
            if (ready > first || first - ready > 4) {
                printf("LDU (DUID %X) with %d frames: ready %d, the first that changes %d\n", duid, n, ready, first);
                bad = 1;
            }
        }
    }
    printf(bad ? "FAILED: p25_tx_ldu_ready()\n" : "p25_tx_ldu_ready() holds for 0 - 9 frames\n");
    return bad;
}

/* The modem, finding the HDU's sync from nothing, gives no air for the 112 dibits from 78 after the sync on (docs of
   OpenDMR, p25.md): given as P25_DIBIT_ERASED, the header is still decoded (13 erasures; its RS(36,20) takes 14,
   two checks kept back); as zeros it is not. 78 .. 221 (16 erasures, after power on): not taken */
static int check_hdu_gap(void)
{
    static uint8_t h[400];
    p25_es es = {{1, 2, 3, 4, 5, 6, 7, 8, 9}, P25_ALGID_CLEAR, 0, 0};
    static const struct {
        int end, mark, flips; /* the air not given: 78 .. end - 1; marked erased; bit errors elsewhere */
        const char *what;
    } cases[] = {{190, 1, 0, "78 .. 189 marked erased"},
                 {190, 0, 0, "78 .. 189 as zeros"},
                 {190, 1, 1, "78 .. 189 marked erased, 1 bit wrong elsewhere"},
                 {222, 1, 0, "78 .. 221 marked erased"},
                 {222, 1, 2, "78 .. 221 marked erased, 2 bits wrong elsewhere"}};
    int bad = 0;
    for (unsigned c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
        int ok = 0, wrong = 0, tries = 200;
        for (int t = 0; t < tries; t++) {
            uint16_t tg = (uint16_t)(rnd() & 0xFFFF);
            int n = p25_tx_hdu(h, NAC, &es, 0, tg, 2);
            for (int f = 0; f < cases[c].flips; f++) { /* in the HDU's codes, outside the gap */
                int at = 60 + (int)(rnd() % 300);
                at = at < 78 ? at : (at - 78 + cases[c].end) % (n - 2);
                h[at] ^= (uint8_t)(1 + (rnd() & 1));
            }
            p25_rx rx;
            p25_rx_init(&rx, P25_NAC_ANY);
            int got = 0;
            for (int i = 0; i < n + 40 && !got; i++) {
                uint8_t d = i < n ? h[i] : (uint8_t)(rnd() & 3);
                if (i >= 78 && i < cases[c].end) {
                    d = cases[c].mark ? (uint8_t)(P25_DIBIT_ERASED | (rnd() & 3)) : 0;
                }
                if (p25_rx_dibit(&rx, d) == P25_DUID_HDU) {
                    got = 1;
                    int right = rx.frame.tgid == tg && memcmp(rx.frame.es.mi, es.mi, 9) == 0;
                    ok += rx.frame.es.ok && right;
                    wrong += rx.frame.es.ok && !right;
                }
            }
        }
        printf("HDU without the air of %s: %d of %d decoded, %d taken wrong\n", cases[c].what, ok, tries, wrong);
        bad |= wrong > 0;
        bad |= (c == 0) && (ok != tries);
    }
    return bad;
}

int main(void)
{
    int n = build();
    int fail = check_ldu_ready();
    fail |= check_hdu_gap();
    printf("a call of %d dibits: HDU, %d x (LDU1, LDU2), TDULC\n", n, PAIRS);
    const double bers[] = {0, 0.005, 0.01, 0.02, 0.05};
    for (unsigned i = 0; i < sizeof(bers) / sizeof(bers[0]); i++) {
        result r = run(n, bers[i]);
        printf("BER %4.1f %%: units %2d (HDU %d, LDU %2d, TDULC %d), LC ok %2d of %d, ES ok %2d of %d, IMBE frames "
               "intact %3d of %d\n",
               bers[i] * 100, r.units, r.hdu, r.ldu, r.tdulc, r.lc_ok, PAIRS + 1, r.es_ok, PAIRS + 1, r.imbe_ok,
               PAIRS * 2 * 9);
        if (bers[i] == 0 && (r.units != 2 + 2 * PAIRS || r.lc_ok != PAIRS + 1 || r.es_ok != PAIRS + 1 ||
                             r.imbe_ok != PAIRS * 2 * 9)) {
            fail = 1;
        }
    }
    printf(fail ? "FAILED: the clean call is not decoded whole\n" : "the clean call is decoded whole\n");
    return fail;
}
