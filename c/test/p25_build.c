/*
 * A P25 call from IMBE frames, for a transmitter that sends whole data units (an MMDVM board, a signal generator):
 * HDU, LDU1 / LDU2 pairs (9 voice frames each, the last LDU padded with silence frames), TDU.
 *
 *   p25_build in.imbe out.bin [nac tgid source] [out.txt]
 *
 * in.imbe: 11-byte frames (the 88 parameter bits u0..u7 MSB first: imbe_ref enc, p25_cli). out.bin: the data units one
 * after the other, each as it goes on the air (frame sync, NID, status dibits included), packed MSB first:
 * HDU 99 bytes, LDU 216 bytes, TDU 18 bytes. out.txt: one line a data unit, its type and its bytes in hex.
 * Defaults: NAC 0x293, group call TG 1 from source 1, clear.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "imbe.h"
#include "p25_frame.h"

#define STATUS 2 /* the status dibits: 10, "unknown, for talkaround" */

static FILE *bin, *txt;
static int units;

static void put(const char *name, const uint8_t *dibits, int n)
{
    uint8_t b[216];
    int bytes = n / 4;
    memset(b, 0, sizeof(b));
    for (int i = 0; i < n; i++) {
        b[i / 4] |= (uint8_t)((dibits[i] & 3) << (6 - 2 * (i % 4)));
    }
    fwrite(b, 1, (size_t)bytes, bin);
    if (txt) {
        fprintf(txt, "%-5s %3d ", name, bytes);
        for (int i = 0; i < bytes; i++) {
            fprintf(txt, "%02X", b[i]);
        }
        fprintf(txt, "\n");
    }
    units++;
}

/* the padding of the last LDU: the input's last frame again (pad the speech with silence before encoding, and that is
   a silence frame) */

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "p25_build in.imbe out.bin [nac tgid source] [out.txt]\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb");
    bin = fopen(argv[2], "wb");
    if (!in || !bin) {
        perror("open");
        return 1;
    }
    uint16_t nac = argc > 3 ? (uint16_t)strtoul(argv[3], NULL, 0) : 0x293;
    uint16_t tg = argc > 4 ? (uint16_t)strtoul(argv[4], NULL, 0) : 1;
    uint32_t src = argc > 5 ? (uint32_t)strtoul(argv[5], NULL, 0) : 1;
    txt = argc > 6 ? fopen(argv[6], "w") : NULL;

    /* the frames */
    uint8_t (*frames)[18] = NULL;
    int nframes = 0, cap = 0;
    uint8_t packed[11];
    while (fread(packed, 1, sizeof(packed), in) == sizeof(packed)) {
        uint8_t u[88];
        for (int i = 0; i < 88; i++) {
            u[i] = (uint8_t)((packed[i >> 3] >> (7 - (i & 7))) & 1);
        }
        if (nframes == cap) {
            cap = cap ? cap * 2 : 256;
            frames = realloc(frames, (size_t)cap * 18);
        }
        imbe_encode_frame(u, frames[nframes++]);
    }
    fclose(in);
    /* whole LDU pairs: 18 frames each, padded with silence */
    int total = ((nframes + 17) / 18) * 18;
    if (total == 0) {
        total = 18;
    }
    frames = realloc(frames, (size_t)total * 18);
    for (int f = nframes; f < total; f++) {
        if (nframes) {
            memcpy(frames[f], frames[nframes - 1], 18);
        } else {
            uint8_t u[88] = {0};
            imbe_encode_frame(u, frames[f]);
        }
    }

    p25_lc lc;
    memset(&lc, 0, sizeof(lc));
    lc.lco = 0x00;
    lc.tgid = tg;
    lc.src = src;
    p25_lc_pack(&lc);
    p25_es es;
    memset(&es, 0, sizeof(es));
    es.algid = P25_ALGID_CLEAR;

    uint8_t d[P25_FRAME_MAX];
    put("HDU", d, p25_tx_hdu(d, nac, &es, 0, tg, STATUS));
    for (int k = 0; k < total / 9; k++) {
        uint8_t lsd[2] = {0, 0};
        int duid = (k & 1) ? P25_DUID_LDU2 : P25_DUID_LDU1;
        put((k & 1) ? "LDU2" : "LDU1", d,
            p25_tx_ldu(d, nac, duid, (const uint8_t(*)[18])frames[k * 9], &lc, &es, lsd, STATUS));
    }
    put("TDU", d, p25_tx_tdu(d, nac, STATUS));
    fclose(bin);
    if (txt) {
        fclose(txt);
    }
    printf("%d voice frames (%.1f s), %d data units: HDU, %d LDUs, TDU; NAC %03X TG %u source %u\n", nframes,
           nframes * 0.02, units, total / 9, nac, tg, (unsigned)src);
    free(frames);
    return 0;
}
