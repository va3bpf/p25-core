/*
 * The error-correcting codes of the P25 Phase 1 common air interface (TIA-102.BAAA), all MSB first: the first bit
 * (or 6-bit symbol) sent is the highest-degree coefficient.
 *
 *   BCH(63,16,23)   NID: NAC (12) + DUID (4), generator 6331141367235453 (octal), roots alpha^1..alpha^22 in
 *                   GF(64) (x^6 + x + 1): corrects 11 bits (Berlekamp-Massey)
 *   Golay(24,12,8)  HDU (shortened to (18,6)) and TDULC words: Golay(23,12) of g = 0xC75 (c/imbe) and an even
 *                   parity bit
 *   Hamming(10,6,3) LDU link control / encryption sync words
 *   Reed-Solomon    over GF(64): (24,12,13) link control, (24,16,9) encryption sync, (36,20,17) header; roots
 *                   alpha^1..alpha^(n-k)
 */
#ifndef P25_FEC_H
#define P25_FEC_H

#include <stdint.h>

/* 16 data bits -> 63-bit codeword (bit 62 first) */
uint64_t p25_bch63_encode(uint16_t data);
/* -> the 16 data bits; returns the bits corrected, -1 if uncorrectable */
int p25_bch63_decode(uint64_t cw, uint16_t *data);

/* 12 data bits -> 24 bits (data, 11 check bits, parity) */
uint32_t p25_golay24_encode(uint32_t d12);
/* -> 12 data bits; returns the bits corrected (up to 3), -1 if 4 errors were detected */
int p25_golay24_decode(uint32_t c24, uint32_t *d12);
/* Golay(18,6): the 6-bit word, then the 12 parity bits of Golay(24,12) on it (its 6 high data bits zero) */
uint32_t p25_golay18_encode(uint32_t d6);
int p25_golay18_decode(uint32_t c18, uint32_t *d6);

/* 6 data bits -> 10 bits (data, 4 parity) */
uint32_t p25_hamming10_encode(uint32_t d6);
/* -> 6 data bits; returns the bits corrected, -1 if the syndrome is no single-bit error */
int p25_hamming10_decode(uint32_t c10, uint32_t *d6);

/* Shortened Reed-Solomon over GF(64): sym[0..n-1] the 6-bit symbols as sent (k data, then n - k parity) */
void p25_rs_encode(uint8_t *sym, int n, int k);
/* corrects sym in place; returns the symbols corrected, -1 if uncorrectable (sym left as it was) */
int p25_rs_decode(uint8_t *sym, int n, int k);
/* The same with erasures: the ne symbols sym[era[0..ne-1]] known to be unreliable (not received). Corrects e errors
   with them while 2 e + ne <= n - k; returns the errors and erasures put right (as located), -1 if uncorrectable. */
int p25_rs_decode_erasures(uint8_t *sym, int n, int k, const uint8_t *era, int ne);

#endif
