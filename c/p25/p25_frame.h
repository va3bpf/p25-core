/*
 * P25 Phase 1 data units (TIA-102.BAAA), receive and transmit, on dibits:
 *
 *   dibit  01 -> +3 (+1800 Hz), 00 -> +1, 10 -> -1, 11 -> -3 (-1800 Hz), the first bit the high one
 *   frame  frame sync (24 dibits), NID (32), the payload; a status dibit after every 35 others (frame dibits 35, 71,
 *          107, ..., counted from the sync)
 *   NID    NAC (12 bits) + DUID (4), BCH(63,16), then one bit set for LDU1 and LDU2 (DSD's table)
 *
 * Lengths with sync and status: HDU 396 dibits, LDU1 / LDU2 864 (180 ms), TDU 72, TDULC 216.
 */
#ifndef P25_FRAME_H
#define P25_FRAME_H

#include <stdint.h>

#define P25_SYNC 0x5575F5FF77FFull

enum {
    P25_DUID_HDU = 0x0,
    P25_DUID_TDU = 0x3,
    P25_DUID_LDU1 = 0x5,
    P25_DUID_TSDU = 0x7,
    P25_DUID_LDU2 = 0xA,
    P25_DUID_PDU = 0xC,
    P25_DUID_TDULC = 0xF,
};

#define P25_ALGID_CLEAR 0x80
#define P25_NAC_ANY 0xF7E   /* the receiver takes any NAC */
#define P25_FRAME_MAX 864   /* dibits, LDU */

/* link control (LDU1, TDULC): the 72 bits, and the fields of the voice formats */
typedef struct {
    uint8_t bytes[9];
    int ok;          /* Reed-Solomon decoded */
    uint8_t lco;     /* link control opcode: 0x00 group voice, 0x03 unit to unit */
    int protect;     /* the LC is encrypted (P bit) */
    uint8_t mfid;
    uint8_t svc;     /* service options: 0x40 encrypted, 0x80 emergency */
    uint16_t tgid;   /* group voice */
    uint32_t dst;    /* unit to unit */
    uint32_t src;
} p25_lc;

/* encryption sync (LDU2, HDU) */
typedef struct {
    uint8_t mi[9];
    uint8_t algid;
    uint16_t kid;
    int ok;
} p25_es;

typedef struct {
    int duid;
    uint16_t nac;
    int nid_errors;
    int fec_errors;       /* bits / symbols the payload's codes corrected */
    uint8_t imbe[9][18];  /* LDU: the 9 voice frames, c0..c7 MSB first (c/imbe's format) */
    p25_lc lc;            /* LDU1, TDULC */
    p25_es es;            /* LDU2, HDU */
    uint8_t mfid;         /* HDU */
    uint16_t tgid;        /* HDU */
    uint8_t lsd[2];       /* LDU low speed data (as received, not corrected) */
    uint8_t status[24];   /* the status dibits */
    int nstatus;
} p25_frame;

typedef struct {
    uint64_t shift;       /* the last 48 bits */
    int sync_tolerance;   /* bits that may differ from the frame sync (4) */
    uint16_t nac;         /* the NAC to take, P25_NAC_ANY for any */
    int in_frame;
    int pos;              /* dibits of the frame so far, sync included */
    int len;              /* the frame's length once the NID is decoded */
    uint8_t d[P25_FRAME_MAX]; /* the frame's dibits without the status dibits, from the NID on (P25_DIBIT_ERASED kept) */
    int nd;
    p25_frame frame;
} p25_rx;

void p25_rx_init(p25_rx *rx, uint16_t nac);
/* or'd into a dibit given to p25_rx_dibit(): not received (the modem gave nothing for it), its value unknown; the
   codes of the data unit take the words it is in as erasures */
#define P25_DIBIT_ERASED 4
/* one received dibit (P25_DIBIT_ERASED: not received); returns the DUID when a data unit is complete (rx->frame),
   else -1 */
int p25_rx_dibit(p25_rx *rx, uint8_t dibit);
/* While an LDU comes in (its NID decoded): how many of its 9 voice frames are complete, from the first; -1: not in an
   LDU. For decoding each as it comes, not 180 ms later with the whole LDU. */
int p25_rx_ldu_frames(const p25_rx *rx);
/* Voice frame k (complete: k < p25_rx_ldu_frames()) of the LDU coming in, as p25_frame imbe[k] */
void p25_rx_ldu_frame(const p25_rx *rx, int k, uint8_t out[18]);

/* ---- transmit: the whole data unit as dibits (sync, NID, status dibits included); returns its length ---- */
int p25_tx_hdu(uint8_t *out, uint16_t nac, const p25_es *es, uint8_t mfid, uint16_t tgid, uint8_t status);
/* LDU1 carries lc->bytes, LDU2 es; lsd may be NULL */
int p25_tx_ldu(uint8_t *out, uint16_t nac, int duid, const uint8_t imbe[9][18], const p25_lc *lc, const p25_es *es,
               const uint8_t lsd[2], uint8_t status);
/* The air dibits of an LDU, from its start, that p25_tx_ldu() gives the same whatever its voice frames from `frames`
   (0 - 9) on: up to where that frame starts (the LC / ES and the low speed data between are known before). For
   sending each part as soon as its voice frame is encoded, not 180 ms later with the whole LDU. */
int p25_tx_ldu_ready(int frames);
int p25_tx_tdu(uint8_t *out, uint16_t nac, uint8_t status);
int p25_tx_tdulc(uint8_t *out, uint16_t nac, const p25_lc *lc, uint8_t status);

/* the voice formats' 72 bits from the fields, and back */
void p25_lc_pack(p25_lc *lc);
void p25_lc_unpack(p25_lc *lc);

#endif
