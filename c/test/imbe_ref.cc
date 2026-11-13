/*
 * Reference IMBE decoder for comparisons on the PC: 11-byte frames (p25_cli's out.imbe) -> WAV, with Pavel Yazev's
 * imbe_vocoder (ref/imbe_vocoder, GPL v3, fetched by fetch_refs.sh). Built only by "c/build.sh ref", for local
 * comparisons: the program links GPL code and is not distributed; nothing of imbe_vocoder goes into this repository.
 *
 *   imbe_ref dec in.imbe out.wav     11-byte frames -> WAV
 *   imbe_ref enc in.raw out.imbe     8 kHz 16-bit mono raw -> 11-byte frames (its encoder)
 */
#include <cstdint>
#include <cstdio>

#include "typedef.h" /* imbe_vocoder's headers: IMBE_PARAM first */
#include "imbe.h"
#include "imbe_vocoder.h"

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
    if (argc < 4) {
        fprintf(stderr, "imbe_ref dec in.imbe out.wav | imbe_ref enc in.raw out.imbe\n");
        return 2;
    }
    FILE *in = fopen(argv[2], "rb");
    FILE *out = fopen(argv[3], "wb");
    if (!in || !out) {
        perror("open");
        return 1;
    }
    imbe_vocoder voc;
    uint8_t frame[11];
    int16_t pcm[160];
    uint32_t samples = 0;
    if (argv[1][0] == 'e') {
        int frames = 0;
        while (fread(pcm, 2, 160, in) == 160) {
            voc.encode_4400(pcm, frame);
            fwrite(frame, 1, sizeof(frame), out);
            frames++;
        }
        fclose(out);
        fclose(in);
        printf("wrote %s (%d frames)\n", argv[3], frames);
        return 0;
    }
    wav_header(out, 0);
    while (fread(frame, 1, sizeof(frame), in) == sizeof(frame)) {
        voc.decode_4400(pcm, frame);
        for (int n = 0; n < 160; n++) {
            put_le(out, (uint16_t)pcm[n], 2);
        }
        samples += 160;
    }
    fseek(out, 0, SEEK_SET);
    wav_header(out, samples);
    fclose(out);
    fclose(in);
    printf("wrote %s (%u frames)\n", argv[3], samples / 160);
    return 0;
}
