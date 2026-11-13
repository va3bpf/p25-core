/*
 * IMBE full-rate (P25 Phase 1, 7200x4400) vocoder in C, fixed point for the radio's CK803S (no FPU).
 *
 * A frame is 20 ms: 88 parameter bits u0..u7 (12 12 12 12 11 11 11 7) and 56 error-correction bits, 144 channel bits:
 * c0..c3 Golay(23,12), c4..c6 Hamming(15,11), c7 as it is, c1..c6 scrambled with a sequence seeded by u0.
 *
 * Frames are 18 bytes, the 144 channel bits MSB first in the order c0 (bit 22 first) .. c7 (bit 6 first), i.e. after
 * the P25 interleave is undone. The parameter bits are one a byte (0 / 1), u0 MSB first .. u7 (mbelib's imbe_d order).
 */
#ifndef IMBE_H
#define IMBE_H

#include <stdint.h>

#define IMBE_N 160 /* samples a frame, 20 ms at 8 kHz */

/* ---- frame layer (imbe_frame.c) ---- */

typedef struct {
    int e0;     /* bit errors corrected in c0 */
    int etotal; /* in c0..c6 */
    int bad;    /* a word with more errors than its code corrects (its decoding is a guess) */
} imbe_errors;

/* 18 bytes -> 88 parameter bits */
void imbe_decode_frame(const uint8_t in18[18], uint8_t u88[88], imbe_errors *err);
/* 88 parameter bits -> 18 bytes */
void imbe_encode_frame(const uint8_t u88[88], uint8_t out18[18]);

/* the codes, for the tests: data in the high bits (Golay 12 of 23, Hamming 11 of 15) */
uint32_t imbe_golay23_encode(uint32_t d12);
uint32_t imbe_golay23_decode(uint32_t c23, int *errors); /* -> 12 data bits; errors -1 if more than 3 */
uint32_t imbe_hamming15_encode(uint32_t d11);
uint32_t imbe_hamming15_decode(uint32_t c15, int *errors); /* -> 11 data bits */

/* ---- parameters (imbe_parms.c) ---- */

typedef struct {
    float w0;
    uint32_t w0q;       /* w0 as a phase step, 2^32 = 2 pi */
    int L;
    int repeat;
    int Vl[58];
    int32_t Ml[58];     /* spectral amplitudes, Q6 */
    int32_t log2Ml[58]; /* Q12 */
    uint32_t PHIl[58];  /* phases, 2^32 = 2 pi */
    uint32_t PSIl[58];
} imbe_parms;

void imbe_parms_init(imbe_parms *p);
/* 0 voice, 1 not a valid frame (b0 above 207: silence or an invalid pitch) */
int imbe_decode_parms(const uint8_t u88[88], imbe_parms *cur, imbe_parms *prev);

/* ---- encoder (imbe_enc.c) ---- */

/* the quantizer: the decoder's last frame, the base of its amplitude prediction */
typedef struct {
    imbe_parms prev;
} imbe_quantizer;

void imbe_quant_init(imbe_quantizer *q);
/* the b0 nearest to a fundamental (2^32 = 2 pi) */
int imbe_b0_of(uint32_t w0q);
/* b0, a voicing decision for each harmonic Vl[1..L] and log2 Ml[1..L] (Q12), L = the L of b0 -> 88 parameter bits */
void imbe_quantize(imbe_quantizer *q, int b0, const uint8_t *Vl, const int32_t *log2Ml, uint8_t u88[88]);

/* the encoder (imbe_ana.c): the analysis of dmr-vocoder's AMBE+2 encoder with IMBE's parameters */
#define IMBE_ENC_HIST 1024 /* input history kept (samples) */
#define IMBE_ENC_DELAY 1   /* frames the output lags the input: push that many frames of zeros at the end */

typedef struct {
    int16_t x[IMBE_ENC_HIST];    /* the latest input samples */
    long n_in;                   /* samples received */
    int32_t dc_x;                /* the input's DC blocker: the last input sample */
    int64_t dc_y;                /* and the last output, Q15 */
    /* pitch tracker */
    int32_t E[3][217];           /* E(P) of frame E_tag[k % 3], Q16 */
    long E_tag[3];
    int hist_i[2];
    int32_t hist_e[2];           /* Q16 */
    int32_t track_e;             /* Q16 */
    float last_w0, prev_p;
    int have_last;
    int prev_bands[16];          /* the last frame's voicing in 500 Hz steps */
    /* voice activity */
    float floor_db, peak_db;
    int have_level, last_raw;
    imbe_quantizer q;
} imbe_encoder;

void imbe_enc_init(imbe_encoder *e);
/* 160 input samples in; 88 parameter bits out (returns 1) once the look-ahead is there */
int imbe_enc_frame(imbe_encoder *e, const int16_t in[IMBE_N], uint8_t u88[88]);
/* the same as an 18-byte frame (with the error-correction, imbe_encode_frame) */
int imbe_enc_frame18(imbe_encoder *e, const int16_t in[IMBE_N], uint8_t out18[18]);

/* ---- decoder (imbe_dec.c) ---- */

typedef struct {
    imbe_parms cur, prev, prev_enh;
    int32_t noise_tail[IMBE_N]; /* Q4 */
    int32_t hp[2][4];           /* the two high-pass sections' x1 x2 y1 y2, Q8 */
    uint32_t rng;
} imbe_decoder;

void imbe_dec_init(imbe_decoder *d);
/* 18 bytes -> 160 samples; returns the frame kind: 'V' voice, 'R' repeat, 'M' muted */
int imbe_dec_frame(imbe_decoder *d, const uint8_t in18[18], int16_t out[IMBE_N]);
/* the same from parameter bits with their error counts (the frame layer's) */
int imbe_dec_parms_frame(imbe_decoder *d, const uint8_t u88[88], const imbe_errors *err, int16_t out[IMBE_N]);

#endif
