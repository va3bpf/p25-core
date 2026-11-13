# P25 Phase 1 core (IMBE vocoder + common air interface)

[中文](README.md) · Used by [OpenDMR](https://github.com/va3bpf/OpenDMR), the open firmware for the Baofeng DM-32UV
(its `p25/` submodule).

The two parts of APCO P25 Phase 1 conventional (simplex / repeater) voice, written from scratch. The C code does not
depend on the radio and is tested on the PC:

- **IMBE vocoder** (full rate 7200x4400: a frame every 20 ms, 88 bits of parameters + 56 bits of error correction =
  144 bits): decoder and encoder, fixed-point C
- **Common air interface**: frame sync, NID (BCH), HDU / LDU1 / LDU2 / TDU / TDULC received and sent, Golay / Hamming /
  Reed-Solomon codes (Reed-Solomon with errors and erasures: symbols not received are taken as erasures), link
  control, encryption sync

The radio's side (raw receive and transmit through the HR_C7000 modem, codeplug, user interface) is in OpenDMR.

**Status**: two DM-32UV talk to each other simplex, and with an MMDVM hotspot, in OpenDMR (see its
[docs/p25.md](https://github.com/va3bpf/OpenDMR/blob/develop/docs/p25.md), in Chinese, with the known issues). The
plan and the protocol parameters checked: [docs/plan.md](docs/plan.md) (Chinese).

## Contents

| Path | What |
|---|---|
| `c/imbe/` | The IMBE vocoder (fixed-point C, decoder and encoder) |
| `c/p25/` | The air interface: `p25_frame.c` frames in and out, `p25_fec.c` the codes |
| `c/test/` | Tests and command line tools for the PC, built by `c/build.sh` (WinLibs MinGW-w64 GCC) |
| `py/` | Python: the IMBE tables, comparisons with reference code, a C4FM demodulator |
| `docs/plan.md` | Plan, references and their licences, protocol parameters checked |
| `ref/` | Reference sources (mbelib, dsd): **not in git**, fetched by `fetch_refs.sh` |
| `bin/`, `out/` | Build and test output, **not in git** |

## Licence

See [LICENSE](LICENSE): free to use, change and share for non-commercial purposes; commercial use needs the author's
written permission. The parts ported from mbelib and DSD are under their ISC licence as well. IMBE is DVSI's
technology: this licence grants no patent rights.
