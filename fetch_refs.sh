#!/bin/sh
# Fetch the reference sources into ref/ (not in git):
#   mbelib  - IMBE 7200x4400 decoder (ISC licence): the parameter decoding, the tables, the Golay / Hamming codes
#   dsd     - P25 Phase 1 frame layout and IMBE interleave schedule (ISC licence; not its GPL D-STAR files, and not
#             its BCH decoder, which is IT++'s, GPL)
set -e
cd "$(dirname "$0")"
mkdir -p ref
cd ref
[ -d mbelib ] || git clone --depth 1 https://github.com/szechyjs/mbelib.git
[ -d dsd ] || git clone --depth 1 https://github.com/szechyjs/dsd.git
# imbe_vocoder - Pavel Yazev's fixed-point IMBE (GPL v3): a reference decoder for comparisons on the PC only
#                (c/build.sh ref); none of it goes into this repository
[ -d imbe_vocoder ] || git clone --depth 1 https://github.com/nostar/imbe_vocoder.git
