#ifndef IMBE_ANA_TABLES_H
#define IMBE_ANA_TABLES_H

/* the hot tables can go into RAM (-DIMBE_RAM_DATA=...) */
#ifndef IMBE_RAM_DATA
#define IMBE_RAM_DATA
#endif

extern const short ImbeACosQ12[1025];      /* cos(2 pi i / 4096) x 4096, a quarter of the circle (cos_q12()) */
extern const short ImbeACos2048Q15[513];   /* cos(2 pi i / 2048) x 32768 (32767 at most), a quarter (cos2048_q15()) */
extern const unsigned short ImbeAExp2Q15[257]; /* 2^(i / 256) x 32768 */
extern const short ImbeALog2Q12[257];      /* log2(1 + i / 256) x 4096 */
extern const short ImbeAHammingQ15[256];   /* the analysis window */
extern const short ImbeAWi2Q15[301];       /* the pitch window w_I^2, / its peak */
extern const short ImbeALpfQ15[21];        /* the 1 kHz low-pass of the pitch */
extern const short ImbeAWinTQ12[2049];     /* the analysis window's transform, real, 32 a bin, / its peak */
extern const int ImbeALogPQ16[217];        /* ln P of the pitch grid P = 20 + i / 2 */
extern const short ImbeAPmQ14[217];        /* P max(w_I^2) */
extern const short ImbeAInvDQ12[217];      /* 1 / (1 - P sum(w_I^4)) */
extern const unsigned char ImbeAPLo[217];  /* the grid indices within -20 % .. +20 % of each P */
extern const unsigned char ImbeAPHi[217];

#endif
