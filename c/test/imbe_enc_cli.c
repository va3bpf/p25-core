/*
 * Our IMBE encoder on the PC: 8 kHz 16-bit mono raw -> 11-byte frames (the 88 parameter bits u0..u7 MSB first, as
 * imbe_ref enc writes them).
 *
 *   imbe_enc_cli in.raw out.imbe
 */
#include <stdio.h>
#include <string.h>

#include "imbe.h"

static imbe_encoder enc;

static void put(FILE *f, const uint8_t u[88])
{
    uint8_t b[11];
    memset(b, 0, sizeof(b));
    for (int i = 0; i < 88; i++) {
        b[i >> 3] |= (uint8_t)(u[i] << (7 - (i & 7)));
    }
    fwrite(b, 1, sizeof(b), f);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "imbe_enc_cli in.raw out.imbe\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb"), *out = fopen(argv[2], "wb");
    if (!in || !out) {
        perror("open");
        return 1;
    }
    imbe_enc_init(&enc);
    int16_t pcm[IMBE_N];
    uint8_t u[88];
    int frames = 0;
    while (fread(pcm, 2, IMBE_N, in) == IMBE_N) {
        if (imbe_enc_frame(&enc, pcm, u)) {
            put(out, u);
            frames++;
        }
    }
    memset(pcm, 0, sizeof(pcm));
    for (int i = 0; i < IMBE_ENC_DELAY; i++) {
        if (imbe_enc_frame(&enc, pcm, u)) {
            put(out, u);
            frames++;
        }
    }
    fclose(in);
    fclose(out);
    printf("wrote %s (%d frames)\n", argv[2], frames);
    return 0;
}
