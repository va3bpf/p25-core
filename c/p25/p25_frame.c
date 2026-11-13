/*
 * P25 Phase 1 data units: receive (sync, NID, status dibits, HDU / LDU1 / LDU2 / TDU / TDULC) and transmit.
 *
 * The positions within the data units follow DSD's p25p1_*.c (ISC licence, see LICENSE), and the IMBE interleave
 * schedule is DSD's p25p1_const.h (iW, iX, iY, iZ). The codes are ours (p25_fec.c).
 *
 * Payload positions, in dibits from the end of the NID with the status dibits taken out:
 *   LDU   IMBE 1 at 0, IMBE 2 at 72, 4 words at 144, IMBE 3 at 164, 4 words at 236, IMBE 4 at 256, 4 words at 328,
 *         IMBE 5 at 348, 4 words at 420, IMBE 6 at 440, 4 words at 512, IMBE 7 at 532, 4 words at 604, IMBE 8 at 624,
 *         low speed data at 696 (byte, its parity, byte, its parity: 16 dibits), IMBE 9 at 712: 784 dibits.
 *         The 24 words (6 bits + Hamming(10,6), 5 dibits) in order are a Reed-Solomon codeword: LDU1 (24,12) the
 *         link control, LDU2 (24,16) the encryption sync.
 *   HDU   36 words of 6 bits + Golay(18,6) parity (9 dibits): RS(36,20), the header; then 5 null dibits
 *   TDULC 12 words of 12 bits + Golay(24,12) parity (12 dibits): RS(24,12), the link control; then 10 null dibits
 *   TDU   14 null dibits
 */
#include <string.h>

#include "p25_fec.h"
#include "p25_frame.h"

#define NID_DIBITS 32

static const uint16_t LDU_IMBE[9] = {0, 72, 164, 256, 348, 440, 532, 624, 712};
static const uint16_t LDU_WORDS[6] = {144, 236, 328, 420, 512, 604};
#define LDU_LSD 696
#define LDU_PAYLOAD 784

/* DSD's IMBE interleave schedule: dibit j carries bit iX[j] of c_iW[j] (its high bit) and bit iZ[j] of c_iY[j] */
static const uint8_t iW[72] = {0, 2, 4, 1, 3, 5, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6,
                               0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 6,
                               0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 4, 7, 0, 3, 5, 2, 4, 7};
static const uint8_t iX[72] = {22, 20, 10, 20, 18, 0,  20, 18, 8,  18, 16, 13, 18, 16, 6,  16, 14, 11,
                               16, 14, 4,  14, 12, 9,  14, 12, 2,  12, 10, 7,  12, 10, 0,  10, 8,  5,
                               10, 8,  13, 8,  6,  3,  8,  6,  11, 6,  4,  1,  6,  4,  9,  4,  2,  6,
                               4,  2,  7,  2,  0,  4,  2,  0,  5,  0,  13, 2,  0,  21, 3,  21, 11, 0};
static const uint8_t iY[72] = {1, 3, 5, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 4,
                               1, 3, 6, 0, 2, 4, 1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5, 1, 3, 6, 0, 2, 5,
                               1, 3, 6, 0, 2, 5, 1, 3, 7, 0, 2, 5, 1, 4, 7, 0, 3, 5, 2, 4, 7, 1, 3, 5};
static const uint8_t iZ[72] = {21, 19, 1,  21, 19, 9,  19, 17, 14, 19, 17, 7,  17, 15, 12, 17, 15, 5,
                               15, 13, 10, 15, 13, 3,  13, 11, 8,  13, 11, 1,  11, 9,  6,  11, 9,  14,
                               9,  7,  4,  9,  7,  12, 7,  5,  2,  7,  5,  10, 5,  3,  0,  5,  3,  8,
                               3,  1,  5,  3,  1,  6,  1,  14, 3,  1,  22, 4,  22, 12, 1,  22, 20, 2};

static const uint8_t CW_LEN[8] = {23, 23, 23, 23, 15, 15, 15, 7};

/* ---- bits ---- */

/* n dibits from d[off], MSB first */
static uint64_t dbits(const uint8_t *d, int off, int n)
{
    uint64_t v = 0;
    for (int i = 0; i < n; i++) {
        v = (v << 2) | (d[off + i] & 3u);
    }
    return v;
}

/* any of n dibits from d[off] not received */
static int erased(const uint8_t *d, int off, int n)
{
    for (int i = 0; i < n; i++) {
        if (d[off + i] & P25_DIBIT_ERASED) {
            return 1;
        }
    }
    return 0;
}

static void put_dbits(uint8_t *d, int off, int n, uint64_t v)
{
    for (int i = 0; i < n; i++) {
        d[off + i] = (uint8_t)((v >> (2 * (n - 1 - i))) & 3);
    }
}

/* 6-bit symbols <-> bytes, MSB first */
static void sym_to_bytes(const uint8_t *sym, int nsym, uint8_t *out)
{
    int nbits = nsym * 6;
    memset(out, 0, (size_t)((nbits + 7) / 8));
    for (int b = 0; b < nbits; b++) {
        if ((sym[b / 6] >> (5 - b % 6)) & 1) {
            out[b >> 3] |= (uint8_t)(0x80 >> (b & 7));
        }
    }
}

static void bytes_to_sym(const uint8_t *in, int nsym, uint8_t *sym)
{
    for (int s = 0; s < nsym; s++) {
        uint8_t v = 0;
        for (int i = 0; i < 6; i++) {
            int b = s * 6 + i;
            v = (uint8_t)((v << 1) | ((in[b >> 3] >> (7 - (b & 7))) & 1));
        }
        sym[s] = v;
    }
}

/* ---- IMBE frames ---- */

static void imbe_deinterleave(const uint8_t *d, uint8_t out18[18])
{
    uint32_t c[8] = {0};
    for (int j = 0; j < 72; j++) {
        c[iW[j]] |= (uint32_t)((d[j] >> 1) & 1) << iX[j];
        c[iY[j]] |= (uint32_t)(d[j] & 1) << iZ[j];
    }
    memset(out18, 0, 18);
    int pos = 0;
    for (int i = 0; i < 8; i++) {
        for (int b = CW_LEN[i] - 1; b >= 0; b--, pos++) {
            if ((c[i] >> b) & 1) {
                out18[pos >> 3] |= (uint8_t)(0x80 >> (pos & 7));
            }
        }
    }
}

static void imbe_interleave(const uint8_t in18[18], uint8_t *d)
{
    uint32_t c[8] = {0};
    int pos = 0;
    for (int i = 0; i < 8; i++) {
        for (int b = CW_LEN[i] - 1; b >= 0; b--, pos++) {
            c[i] |= (uint32_t)((in18[pos >> 3] >> (7 - (pos & 7))) & 1) << b;
        }
    }
    for (int j = 0; j < 72; j++) {
        d[j] = (uint8_t)((((c[iW[j]] >> iX[j]) & 1) << 1) | ((c[iY[j]] >> iZ[j]) & 1));
    }
}

/* ---- link control ---- */

void p25_lc_unpack(p25_lc *lc)
{
    const uint8_t *b = lc->bytes;
    lc->protect = (b[0] >> 7) & 1;
    lc->lco = b[0] & 0x3F;
    lc->mfid = b[1];
    lc->svc = b[2];
    lc->tgid = 0;
    lc->dst = 0;
    if (lc->lco == 0x03) {
        lc->dst = ((uint32_t)b[3] << 16) | ((uint32_t)b[4] << 8) | b[5];
    } else {
        lc->tgid = (uint16_t)((b[4] << 8) | b[5]);
    }
    lc->src = ((uint32_t)b[6] << 16) | ((uint32_t)b[7] << 8) | b[8];
}

void p25_lc_pack(p25_lc *lc)
{
    uint8_t *b = lc->bytes;
    memset(b, 0, 9);
    b[0] = (uint8_t)((lc->protect ? 0x80 : 0) | (lc->lco & 0x3F));
    b[1] = lc->mfid;
    b[2] = lc->svc;
    if (lc->lco == 0x03) {
        b[3] = (uint8_t)(lc->dst >> 16);
        b[4] = (uint8_t)(lc->dst >> 8);
        b[5] = (uint8_t)lc->dst;
    } else {
        b[4] = (uint8_t)(lc->tgid >> 8);
        b[5] = (uint8_t)lc->tgid;
    }
    b[6] = (uint8_t)(lc->src >> 16);
    b[7] = (uint8_t)(lc->src >> 8);
    b[8] = (uint8_t)lc->src;
}

static void es_from_bytes(p25_es *es, const uint8_t *b)
{
    memcpy(es->mi, b, 9);
    es->algid = b[9];
    es->kid = (uint16_t)((b[10] << 8) | b[11]);
}

static void es_to_bytes(const p25_es *es, uint8_t *b)
{
    memcpy(b, es->mi, 9);
    b[9] = es->algid;
    b[10] = (uint8_t)(es->kid >> 8);
    b[11] = (uint8_t)es->kid;
}

/* ---- receive ---- */

void p25_rx_init(p25_rx *rx, uint16_t nac)
{
    memset(rx, 0, sizeof(*rx));
    rx->sync_tolerance = 4;
    rx->nac = nac;
}

static int frame_length(int duid)
{
    switch (duid) {
    case P25_DUID_HDU:
        return 396;
    case P25_DUID_LDU1:
    case P25_DUID_LDU2:
        return 864;
    case P25_DUID_TDU:
        return 72;
    case P25_DUID_TDULC:
        return 216;
    default:
        return 0; /* TSDU, PDU: variable, not handled */
    }
}

/* Reed-Solomon with the erasures of the words given: two of the n - k checks kept back to find a wrong correction
   (with them all used, any error besides gives a wrong word taken as right), so 2 e + erasures <= n - k - 2 */
static int rs_decode(uint8_t *sym, int n, int k, const uint8_t *era, int ne)
{
    if (ne == 0) {
        return p25_rs_decode(sym, n, k);
    }
    int r = p25_rs_decode_erasures(sym, n, k, era, ne);
    if (r < 0 || 2 * (r - ne) + ne > n - k - 2) {
        return -1;
    }
    return r;
}

static void parse_ldu(p25_rx *rx)
{
    p25_frame *f = &rx->frame;
    const uint8_t *p = rx->d + NID_DIBITS;
    uint8_t sym[24], bytes[12], era[24];
    int ne = 0;
    for (int k = 0; k < 9; k++) {
        imbe_deinterleave(p + LDU_IMBE[k], f->imbe[k]);
    }
    for (int t = 0; t < 24; t++) {
        uint32_t v;
        int at = LDU_WORDS[t / 4] + (t % 4) * 5;
        int e = p25_hamming10_decode((uint32_t)dbits(p, at, 5), &v);
        sym[t] = (uint8_t)v;
        if (erased(p, at, 5)) {
            era[ne++] = (uint8_t)t;
        } else {
            f->fec_errors += e > 0 ? e : 0;
        }
    }
    int k = f->duid == P25_DUID_LDU1 ? 12 : 16;
    int r = rs_decode(sym, 24, k, era, ne);
    f->fec_errors += r > 0 ? r : 0;
    sym_to_bytes(sym, k, bytes);
    if (f->duid == P25_DUID_LDU1) {
        memcpy(f->lc.bytes, bytes, 9);
        f->lc.ok = r >= 0;
        p25_lc_unpack(&f->lc);
    } else {
        es_from_bytes(&f->es, bytes);
        f->es.ok = r >= 0;
    }
    f->lsd[0] = (uint8_t)dbits(p, LDU_LSD, 4);
    f->lsd[1] = (uint8_t)dbits(p, LDU_LSD + 8, 4);
}

static void parse_hdu(p25_rx *rx)
{
    p25_frame *f = &rx->frame;
    const uint8_t *p = rx->d + NID_DIBITS;
    uint8_t sym[36], b[15], era[36];
    int ne = 0;
    for (int i = 0; i < 36; i++) {
        uint32_t v;
        int e = p25_golay18_decode((uint32_t)dbits(p, i * 9, 9), &v);
        sym[i] = (uint8_t)v;
        if (erased(p, i * 9, 9)) {
            era[ne++] = (uint8_t)i;
        } else {
            f->fec_errors += e > 0 ? e : 0;
        }
    }
    int r = rs_decode(sym, 36, 20, era, ne);
    f->fec_errors += r > 0 ? r : 0;
    sym_to_bytes(sym, 20, b);
    memcpy(f->es.mi, b, 9);
    f->mfid = b[9];
    f->es.algid = b[10];
    f->es.kid = (uint16_t)((b[11] << 8) | b[12]);
    f->tgid = (uint16_t)((b[13] << 8) | b[14]);
    f->es.ok = r >= 0;
}

static void parse_tdulc(p25_rx *rx)
{
    p25_frame *f = &rx->frame;
    const uint8_t *p = rx->d + NID_DIBITS;
    uint8_t sym[24], era[24];
    int ne = 0;
    for (int i = 0; i < 12; i++) {
        uint32_t v;
        int e = p25_golay24_decode((uint32_t)dbits(p, i * 12, 12), &v);
        sym[2 * i] = (uint8_t)(v >> 6);
        sym[2 * i + 1] = (uint8_t)(v & 63);
        if (erased(p, i * 12, 12)) {
            era[ne++] = (uint8_t)(2 * i);
            era[ne++] = (uint8_t)(2 * i + 1);
        } else {
            f->fec_errors += e > 0 ? e : 0;
        }
    }
    int r = rs_decode(sym, 24, 12, era, ne);
    f->fec_errors += r > 0 ? r : 0;
    sym_to_bytes(sym, 12, f->lc.bytes);
    f->lc.ok = r >= 0;
    p25_lc_unpack(&f->lc);
}

static void hunt(p25_rx *rx)
{
    rx->in_frame = 0;
    rx->shift = 0;
}

int p25_rx_ldu_frames(const p25_rx *rx)
{
    if (!rx->in_frame || !rx->len ||
        ((rx->frame.duid != P25_DUID_LDU1) && (rx->frame.duid != P25_DUID_LDU2))) {
        return -1;
    }
    int n = 0;
    while ((n < 9) && (rx->nd >= NID_DIBITS + LDU_IMBE[n] + 72)) {
        n++;
    }
    return n;
}

void p25_rx_ldu_frame(const p25_rx *rx, int k, uint8_t out[18])
{
    imbe_deinterleave(rx->d + NID_DIBITS + LDU_IMBE[k], out);
}

int p25_rx_dibit(p25_rx *rx, uint8_t dibit)
{
    uint8_t lost = dibit & P25_DIBIT_ERASED;
    dibit &= 3;
    if (!rx->in_frame) {
        rx->shift = ((rx->shift << 2) | dibit) & 0xFFFFFFFFFFFFull;
        if (__builtin_popcountll(rx->shift ^ P25_SYNC) <= rx->sync_tolerance) {
            rx->in_frame = 1;
            rx->pos = 24;
            rx->nd = 0;
            rx->len = 0;
            memset(&rx->frame, 0, sizeof(rx->frame));
        }
        return -1;
    }
    p25_frame *f = &rx->frame;
    if (rx->pos % 36 == 35) {
        if (f->nstatus < 24) {
            f->status[f->nstatus++] = dibit;
        }
    } else if (rx->nd < P25_FRAME_MAX) {
        rx->d[rx->nd++] = dibit | lost;
    }
    rx->pos++;
    if (!rx->len && rx->nd == NID_DIBITS) {
        uint64_t nid = dbits(rx->d, 0, NID_DIBITS);
        uint16_t v;
        int e = p25_bch63_decode(nid >> 1, &v);
        if (e < 0 || (rx->nac != P25_NAC_ANY && (v >> 4) != rx->nac)) {
            hunt(rx);
            return -1;
        }
        f->nac = v >> 4;
        f->duid = v & 15;
        f->nid_errors = e;
        rx->len = frame_length(f->duid);
        if (!rx->len) {
            hunt(rx); /* reported, but its payload is not read */
            return f->duid;
        }
    }
    if (rx->len && rx->pos >= rx->len) {
        switch (f->duid) {
        case P25_DUID_LDU1:
        case P25_DUID_LDU2:
            parse_ldu(rx);
            break;
        case P25_DUID_HDU:
            parse_hdu(rx);
            break;
        case P25_DUID_TDULC:
            parse_tdulc(rx);
            break;
        default:
            break;
        }
        hunt(rx);
        return f->duid;
    }
    return -1;
}

/* ---- transmit ---- */

/* sync + NID into d[0..55]; returns 56 */
static int head(uint8_t *d, uint16_t nac, int duid)
{
    put_dbits(d, 0, 24, P25_SYNC);
    uint64_t cw = p25_bch63_encode((uint16_t)(((nac & 0xFFF) << 4) | (duid & 15)));
    uint64_t nid = (cw << 1) | (uint64_t)(duid == P25_DUID_LDU1 || duid == P25_DUID_LDU2);
    put_dbits(d, 24, NID_DIBITS, nid);
    return 24 + NID_DIBITS;
}

/* the data dibits u[0..nu-1] (sync on, without status) into out with a status dibit after every 35; padded with
   zeros to len */
static int with_status(uint8_t *out, const uint8_t *u, int nu, int len, uint8_t status)
{
    int j = 0;
    for (int i = 0; i < len; i++) {
        out[i] = (i % 36 == 35) ? (uint8_t)(status & 3) : (j < nu ? u[j++] : 0);
    }
    return len;
}

int p25_tx_ldu(uint8_t *out, uint16_t nac, int duid, const uint8_t imbe[9][18], const p25_lc *lc, const p25_es *es,
               const uint8_t lsd[2], uint8_t status)
{
    uint8_t u[P25_FRAME_MAX], sym[24], bytes[12];
    int n = head(u, nac, duid);
    uint8_t *p = u + n;
    memset(p, 0, LDU_PAYLOAD);
    for (int k = 0; k < 9; k++) {
        imbe_interleave(imbe[k], p + LDU_IMBE[k]);
    }
    int k = duid == P25_DUID_LDU1 ? 12 : 16;
    if (duid == P25_DUID_LDU1) {
        memcpy(bytes, lc->bytes, 9);
    } else {
        es_to_bytes(es, bytes);
    }
    bytes_to_sym(bytes, k, sym);
    p25_rs_encode(sym, 24, k);
    for (int t = 0; t < 24; t++) {
        put_dbits(p, LDU_WORDS[t / 4] + (t % 4) * 5, 5, p25_hamming10_encode(sym[t]));
    }
    /* low speed data; its cyclic (16,8) parity not done yet (zero) */
    put_dbits(p, LDU_LSD, 4, lsd ? lsd[0] : 0);
    put_dbits(p, LDU_LSD + 8, 4, lsd ? lsd[1] : 0);
    return with_status(out, u, n + LDU_PAYLOAD, 864, status);
}

int p25_tx_ldu_ready(int frames)
{
    if (frames >= 9) {
        return 864;
    }
    int e = 24 + NID_DIBITS + LDU_IMBE[frames < 0 ? 0 : frames]; /* data dibits, sync on */
    return e + e / 35;                                            /* with the status dibits among them */
}

int p25_tx_hdu(uint8_t *out, uint16_t nac, const p25_es *es, uint8_t mfid, uint16_t tgid, uint8_t status)
{
    uint8_t u[400], sym[36], b[15];
    int n = head(u, nac, P25_DUID_HDU);
    memcpy(b, es->mi, 9);
    b[9] = mfid;
    b[10] = es->algid;
    b[11] = (uint8_t)(es->kid >> 8);
    b[12] = (uint8_t)es->kid;
    b[13] = (uint8_t)(tgid >> 8);
    b[14] = (uint8_t)tgid;
    bytes_to_sym(b, 20, sym);
    p25_rs_encode(sym, 36, 20);
    for (int i = 0; i < 36; i++) {
        put_dbits(u, n + i * 9, 9, p25_golay18_encode(sym[i]));
    }
    return with_status(out, u, n + 324, 396, status);
}

int p25_tx_tdulc(uint8_t *out, uint16_t nac, const p25_lc *lc, uint8_t status)
{
    uint8_t u[220], sym[24];
    int n = head(u, nac, P25_DUID_TDULC);
    bytes_to_sym(lc->bytes, 12, sym);
    p25_rs_encode(sym, 24, 12);
    for (int i = 0; i < 12; i++) {
        put_dbits(u, n + i * 12, 12, p25_golay24_encode(((uint32_t)sym[2 * i] << 6) | sym[2 * i + 1]));
    }
    return with_status(out, u, n + 144, 216, status);
}

int p25_tx_tdu(uint8_t *out, uint16_t nac, uint8_t status)
{
    uint8_t u[72];
    int n = head(u, nac, P25_DUID_TDU);
    return with_status(out, u, n, 72, status);
}
