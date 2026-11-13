/*
 * Our IMBE decoder on the PC: 11-byte frames (the 88 parameter bits u0 .. u7 MSB first) -> WAV (8 kHz).
 *
 *   imbe_cli in.imbe out.wav
 */
#include <stdio.h>

#include "imbe.h"

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
    if (argc < 3) {
        fprintf(stderr, "imbe_cli in.imbe out.wav\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb"), *out = fopen(argv[2], "wb");
    if (!in || !out) {
        perror("open");
        return 1;
    }
    static imbe_decoder dec;
    imbe_dec_init(&dec);
    uint8_t packed[11], u[88];
    int16_t pcm[IMBE_N];
    uint32_t samples = 0;
    imbe_errors none = {0, 0, 0};
    wav_header(out, 0);
    while (fread(packed, 1, sizeof(packed), in) == sizeof(packed)) {
        for (int i = 0; i < 88; i++) {
            u[i] = (uint8_t)((packed[i >> 3] >> (7 - (i & 7))) & 1);
        }
        imbe_dec_parms_frame(&dec, u, &none, pcm);
        for (int n = 0; n < IMBE_N; n++) {
            put_le(out, (uint16_t)pcm[n], 2);
        }
        samples += IMBE_N;
    }
    fseek(out, 0, SEEK_SET);
    wav_header(out, samples);
    fclose(out);
    fclose(in);
    printf("wrote %s (%u frames)\n", argv[2], samples / IMBE_N);
    return 0;
}
