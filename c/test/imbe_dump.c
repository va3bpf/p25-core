/*
 * The parameters of each frame of a file (11-byte frames): b0, L, voiced harmonics, mean log2 Ml of the voiced and of
 * the unvoiced ones (Q12 -> float), one line a frame. For the encoder's tests (py/enc_test.py).
 *
 *   imbe_dump in.imbe
 */
#include <stdio.h>

#include "imbe.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "imbe_dump in.imbe\n");
        return 2;
    }
    FILE *in = fopen(argv[1], "rb");
    if (!in) {
        perror(argv[1]);
        return 1;
    }
    imbe_parms cur, prev;
    imbe_parms_init(&prev);
    uint8_t packed[11], u[88];
    while (fread(packed, 1, 11, in) == 11) {
        for (int i = 0; i < 88; i++) {
            u[i] = (uint8_t)((packed[i >> 3] >> (7 - (i & 7))) & 1);
        }
        int b0 = (u[0] << 7) | (u[1] << 6) | (u[2] << 5) | (u[3] << 4) | (u[4] << 3) | (u[5] << 2) | (u[85] << 1) | u[86];
        cur = prev;
        if (imbe_decode_parms(u, &cur, &prev)) {
            printf("%d 0 0 0 0\n", b0);
            continue;
        }
        prev = cur;
        int nv = 0, nu = 0;
        double sv = 0, su = 0;
        for (int l = 1; l <= cur.L; l++) {
            if (cur.Vl[l]) {
                nv++;
                sv += cur.log2Ml[l] / 4096.0;
            } else {
                nu++;
                su += cur.log2Ml[l] / 4096.0;
            }
        }
        printf("%d %d %d %.3f %.3f\n", b0, cur.L, nv, nv ? sv / nv : 0, nu ? su / nu : 0);
    }
    return 0;
}
