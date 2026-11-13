/*
 * The quantizer against the decoder: each frame of a file (11-byte frames, u0..u7 MSB first) decoded to its
 * parameters and quantized again must give the same 88 bits (its values are the decoder's reconstruction values).
 *
 *   imbe_requant in.imbe
 */
#include <stdio.h>
#include <string.h>

#include "imbe.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "imbe_requant in.imbe\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb");
    if (!in) {
        perror(argv[1]);
        return 1;
    }
    imbe_parms cur, prev;
    imbe_parms_init(&prev);
    imbe_quantizer q;
    imbe_quant_init(&q);
    uint8_t packed[11];
    int frames = 0, same = 0, skipped = 0, diff_bits = 0, region[8] = {0};
    static const int start[9] = {0, 12, 24, 36, 48, 59, 70, 81, 88};
    while (fread(packed, 1, 11, in) == 11) {
        uint8_t u[88], v[88];
        for (int i = 0; i < 88; i++) {
            u[i] = (uint8_t)((packed[i >> 3] >> (7 - (i & 7))) & 1);
        }
        cur = prev;
        if (imbe_decode_parms(u, &cur, &prev)) {
            skipped++;
            continue;
        }
        prev = cur;
        int b0 = (u[0] << 7) | (u[1] << 6) | (u[2] << 5) | (u[3] << 4) | (u[4] << 3) | (u[5] << 2) | (u[85] << 1) | u[86];
        uint8_t Vl[58];
        for (int l = 0; l < 58; l++) {
            Vl[l] = (uint8_t)cur.Vl[l];
        }
        imbe_quantize(&q, b0, Vl, cur.log2Ml, v);
        int d = 0;
        for (int i = 0; i < 87; i++) { /* u7's last bit is not a parameter */
            if (u[i] != v[i]) {
                d++;
                for (int r = 0; r < 8; r++) {
                    if (i >= start[r] && i < start[r + 1]) {
                        region[r]++;
                    }
                }
            }
        }
        frames++;
        same += d == 0;
        diff_bits += d;
    }
    printf("%d frames (%d not voice frames skipped): %d the same, %d bits differ (u0..u7: %d %d %d %d %d %d %d %d)\n",
           frames, skipped, same, diff_bits, region[0], region[1], region[2], region[3], region[4], region[5],
           region[6], region[7]);
    return 0;
}
