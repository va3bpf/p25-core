/*
 * P25 receive on the PC: dibits (py/c4fm_demod.py) -> a line a data unit, and the voice as a WAV (8 kHz).
 *
 *   p25_cli in.dibits [out.wav [out.imbe]]
 *
 * out.imbe: the voice frames' 88 parameter bits after the error correction, 11 bytes a frame (u0 .. u7 MSB first, the
 * input of other IMBE decoders, e.g. imbe_vocoder's decode_4400).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "imbe.h"
#include "p25_frame.h"

static const char *duid_name(int d)
{
    switch (d) {
    case P25_DUID_HDU:
        return "HDU  ";
    case P25_DUID_TDU:
        return "TDU  ";
    case P25_DUID_LDU1:
        return "LDU1 ";
    case P25_DUID_TSDU:
        return "TSDU ";
    case P25_DUID_LDU2:
        return "LDU2 ";
    case P25_DUID_PDU:
        return "PDU  ";
    case P25_DUID_TDULC:
        return "TDULC";
    default:
        return "?    ";
    }
}

static void put_le(FILE *f, uint32_t v, int n)
{
    for (int i = 0; i < n; i++) {
        fputc((int)((v >> (8 * i)) & 0xFF), f);
    }
}

static void wav_header(FILE *f, uint32_t samples)
{
    fwrite("RIFF", 1, 4, f);
    put_le(f, 36 + samples * 2, 4);
    fwrite("WAVEfmt ", 1, 8, f);
    put_le(f, 16, 4);
    put_le(f, 1, 2);
    put_le(f, 1, 2);
    put_le(f, 8000, 4);
    put_le(f, 16000, 4);
    put_le(f, 2, 2);
    put_le(f, 16, 2);
    fwrite("data", 1, 4, f);
    put_le(f, samples * 2, 4);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "p25_cli in.dibits [out.wav]\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb");
    if (!in) {
        perror(argv[1]);
        return 1;
    }
    FILE *wav = argc > 2 ? fopen(argv[2], "wb") : NULL;
    FILE *frames = argc > 3 ? fopen(argv[3], "wb") : NULL;
    uint32_t samples = 0;
    if (wav) {
        wav_header(wav, 0);
    }
    static p25_rx rx;
    static imbe_decoder dec;
    p25_rx_init(&rx, P25_NAC_ANY);
    imbe_dec_init(&dec);
    long pos = 0;
    int c, counts[16] = {0}, nid_err_total = 0, units = 0;
    int kinds[3] = {0}; /* V R M */
    while ((c = fgetc(in)) != EOF) {
        pos++;
        int duid = p25_rx_dibit(&rx, (uint8_t)c);
        if (duid < 0) {
            continue;
        }
        const p25_frame *f = &rx.frame;
        units++;
        counts[duid & 15]++;
        nid_err_total += f->nid_errors;
        printf("%7.3f s  %s NAC %03X  NID errors %2d", pos / 4800.0, duid_name(duid), f->nac, f->nid_errors);
        if (duid == P25_DUID_LDU1 || duid == P25_DUID_TDULC) {
            printf("  LC %s lco %02X mfid %02X svc %02X tg %u src %u dst %u%s", f->lc.ok ? "ok " : "BAD", f->lc.lco,
                   f->lc.mfid, f->lc.svc, f->lc.tgid, (unsigned)f->lc.src, (unsigned)f->lc.dst,
                   f->lc.protect ? " protected" : "");
        }
        if (duid == P25_DUID_LDU2 || duid == P25_DUID_HDU) {
            printf("  ES %s algid %02X kid %04X", f->es.ok ? "ok " : "BAD", f->es.algid, f->es.kid);
        }
        if (duid == P25_DUID_HDU) {
            printf(" mfid %02X tg %u", f->mfid, f->tgid);
        }
        if (duid == P25_DUID_LDU1 || duid == P25_DUID_LDU2 || duid == P25_DUID_HDU || duid == P25_DUID_TDULC) {
            printf("  fec %d", f->fec_errors);
        }
        if (duid == P25_DUID_LDU1 || duid == P25_DUID_LDU2) {
            printf("  imbe");
            for (int k = 0; k < 9; k++) {
                int16_t pcm[IMBE_N];
                uint8_t u[88];
                imbe_errors e;
                imbe_decode_frame(f->imbe[k], u, &e);
                if (frames) {
                    uint8_t packed[11] = {0};
                    for (int i = 0; i < 88; i++) {
                        packed[i >> 3] |= (uint8_t)(u[i] << (7 - (i & 7)));
                    }
                    fwrite(packed, 1, sizeof(packed), frames);
                }
                int kind = imbe_dec_parms_frame(&dec, u, &e, pcm);
                kinds[kind == 'V' ? 0 : kind == 'R' ? 1 : 2]++;
                printf(" %d%c", e.etotal, kind == 'V' ? ' ' : kind);
                if (wav) {
                    for (int n = 0; n < IMBE_N; n++) {
                        put_le(wav, (uint16_t)pcm[n], 2);
                    }
                    samples += IMBE_N;
                }
            }
        }
        if (duid == P25_DUID_HDU || duid == P25_DUID_TDU || duid == P25_DUID_TDULC) {
            imbe_dec_init(&dec);
        }
        printf("\n");
    }
    printf("\n%d data units:", units);
    for (int d = 0; d < 16; d++) {
        if (counts[d]) {
            printf(" %s %d", duid_name(d), counts[d]);
        }
    }
    printf("; NID bits corrected %d\n", nid_err_total);
    if (kinds[0] + kinds[1] + kinds[2]) {
        printf("IMBE frames: %d voice, %d repeated, %d muted\n", kinds[0], kinds[1], kinds[2]);
    }
    if (wav) {
        fseek(wav, 0, SEEK_SET);
        wav_header(wav, samples);
        fclose(wav);
        printf("wrote %s (%.1f s)\n", argv[2], samples / 8000.0);
    }
    if (frames) {
        fclose(frames);
    }
    fclose(in);
    return 0;
}
