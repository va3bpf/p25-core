#!/bin/sh
# Build the PC tools in bin/ with the GCC of WinLibs MinGW-w64 (toolchain/winlibs of OpenDMR, or GCC on the PATH):
#   imbe_test.exe   the tests of c/imbe (our code only)
#   p25_test.exe    the tests of c/p25
# c/build.sh mbelib: imbe_test.exe with the comparisons against mbelib too (ref/mbelib, fetched by fetch_refs.sh,
# compiled into the test only)
set -e
cd "$(dirname "$0")/.."
for d in ../toolchain/winlibs/mingw64/bin toolchain/winlibs/mingw64/bin; do
    [ -d "$d" ] && export PATH="$(cd "$d" && pwd):$PATH"
done
mkdir -p bin out/obj
IMBE="c/imbe/imbe_frame.c c/imbe/imbe_parms.c c/imbe/imbe_dec.c c/imbe/imbe_enc.c c/imbe/imbe_ana.c c/imbe/imbe_tables.c c/imbe/imbe_ana_tables.c"
CFLAGS="-O2 -Wall -Wextra -std=c11"
if [ "$1" = "ref" ]; then
    # imbe_ref.exe: Pavel Yazev's imbe_vocoder (ref/imbe_vocoder, GPL v3) as a reference decoder, local use only
    R=ref/imbe_vocoder
    for f in $R/*.cc; do
        g++ -std=c++11 -O2 -w -I$R -c -o out/obj/ref_$(basename $f .cc).o $f
    done
    g++ -std=c++11 -O2 -Wall -static -I$R -o bin/imbe_ref.exe c/test/imbe_ref.cc out/obj/ref_*.o
    echo "built bin/imbe_ref.exe"
    exit 0
fi
if [ "$1" = "mbelib" ]; then
    M=ref/mbelib
    for f in mbelib ecc imbe7200x4400 imbe7100x4400 ambe3600x2400 ambe3600x2450; do
        gcc -O2 -w -I$M -c -o out/obj/mbe_$f.o $M/$f.c
    done
    gcc $CFLAGS -DWITH_MBELIB -static -Ic/imbe -I$M -o bin/imbe_test.exe c/test/imbe_test.c $IMBE out/obj/mbe_*.o -lm
else
    gcc $CFLAGS -static -Ic/imbe -o bin/imbe_test.exe c/test/imbe_test.c $IMBE -lm
fi
echo "built bin/imbe_test.exe"
P25="c/p25/p25_fec.c c/p25/p25_frame.c"
gcc $CFLAGS -static -Ic/imbe -Ic/p25 -o bin/p25_fec_test.exe c/test/p25_fec_test.c $P25 $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -Ic/p25 -o bin/p25_frame_test.exe c/test/p25_frame_test.c $P25 $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -Ic/p25 -o bin/p25_cli.exe c/test/p25_cli.c $P25 $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -Ic/p25 -o bin/p25_build.exe c/test/p25_build.c $P25 $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -o bin/imbe_cli.exe c/test/imbe_cli.c $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -o bin/imbe_requant.exe c/test/imbe_requant.c $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -o bin/imbe_enc_cli.exe c/test/imbe_enc_cli.c $IMBE -lm
gcc $CFLAGS -static -Ic/imbe -o bin/imbe_dump.exe c/test/imbe_dump.c $IMBE -lm
echo "built bin/p25_fec_test.exe bin/p25_frame_test.exe bin/p25_cli.exe"
