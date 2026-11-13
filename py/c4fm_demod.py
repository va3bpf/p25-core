"""C4FM (P25 Phase 1) demodulator for recordings on the PC: I/Q WAV (or a discriminator WAV) -> dibits.

    py py/c4fm_demod.py in.wav out.dibits [--disc] [--invert]

I/Q: 2 channels (16 or 24 bit), any rate that is a whole multiple of 4800 * 2 or so (24 kHz: 5 samples a
symbol). Steps: channel low-pass, FM discriminator, carrier offset removed, integrate-and-dump over one symbol (the
C4FM receive filter), Gardner symbol timing, a 4-level slicer that follows the levels and the DC.
--disc: the WAV is already a discriminator output (mono). --invert: flip the polarity (some receivers do).

out.dibits: one dibit a byte (0..3: 01 +3, 00 +1, 10 -1, 11 -3), for c/test/p25_cli.c. Prints how many frame syncs
each polarity shows.
"""

import argparse
import struct

import numpy as np

SYMBOL_RATE = 4800
SYNC = 0x5575F5FF77FF
LEVEL_TO_DIBIT = {3: 1, 1: 0, -1: 2, -3: 3}


def read_wav(path):
    b = open(path, "rb").read()
    if b[:4] != b"RIFF" or b[8:12] != b"WAVE":
        raise SystemExit("not a WAV file")
    pos, fmt, data = 12, None, None
    while pos + 8 <= len(b):
        cid, size = b[pos:pos + 4], struct.unpack("<I", b[pos + 4:pos + 8])[0]
        if cid == b"fmt ":
            fmt = struct.unpack("<HHIIHH", b[pos + 8:pos + 24])
        elif cid == b"data":
            data = b[pos + 8:pos + 8 + size]
            break
        pos += 8 + size + (size & 1)
    if fmt is None or data is None:
        raise SystemExit("no fmt / data chunk")
    tag, ch, fs, _, _, bits = fmt
    if bits == 24:
        a = np.frombuffer(data[:len(data) // 3 * 3], dtype=np.uint8).reshape(-1, 3).astype(np.int32)
        v = a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)
        v = np.where(v & 0x800000, v - (1 << 24), v).astype(np.float64)
    elif bits == 16:
        v = np.frombuffer(data[:len(data) // 2 * 2], dtype="<i2").astype(np.float64)
    else:
        raise SystemExit("%d-bit samples not supported" % bits)
    v = v[:len(v) // ch * ch].reshape(-1, ch)
    return fs, v


def lowpass(x, fs, cutoff, taps=63):
    n = np.arange(taps) - (taps - 1) / 2
    h = np.sinc(2 * cutoff / fs * n) * np.hamming(taps)
    return np.convolve(x, h / h.sum(), mode="same")


def discriminate(iq, fs):
    x = lowpass(iq, fs, 6250)
    d = np.angle(x[1:] * np.conj(x[:-1])) * fs / (2 * np.pi)
    amp = np.abs(x[1:])
    # the carrier offset: the mean frequency where the signal is strong (C4FM's levels average out)
    strong = amp > 0.3 * np.percentile(amp, 90)
    offset = np.mean(d[strong]) if strong.any() else 0.0
    d = d - offset
    d[~strong] = 0.0  # no signal: no frequency
    return d, offset, strong.mean()


def symbols(d, fs):
    """Integrate-and-dump, Gardner timing: one soft value a symbol (Hz)."""
    sps = fs / SYMBOL_RATE
    k = int(round(sps))
    y = np.convolve(d, np.ones(k) / k, mode="same")
    out = []
    t = sps
    gain = 0.02
    prev = 0.0
    while t + sps + 1 < len(y):
        i = int(t)
        f = t - i
        cur = y[i] * (1 - f) + y[i + 1] * f
        tm = t - sps / 2
        im = int(tm)
        fm = tm - im
        mid = y[im] * (1 - fm) + y[im + 1] * fm
        err = (cur - prev) * mid
        scale = np.mean(np.abs(out[-200:])) if out else 0.0
        step = 0.0
        if scale > 1.0:  # no signal (the discriminator zeroed): hold the timing
            step = max(-0.3, min(0.3, gain * err / (scale * scale) * sps))
        t += sps - step
        out.append(cur)
        prev = cur
    return np.array(out)


def slice_levels(s):
    """4-level decisions with the outer level and the DC followed; returns the levels (-3, -1, 1, 3)."""
    unit = np.percentile(np.abs(s[s != 0]), 75) / 3 if np.any(s != 0) else 1.0
    dc = 0.0
    lev = np.zeros(len(s), dtype=np.int8)
    for i, v in enumerate(s):
        x = (v - dc) / unit
        l = 3 if x >= 2 else 1 if x >= 0 else -1 if x >= -2 else -3
        lev[i] = l
        if v != 0:
            dc += 0.002 * (x - l) * unit
            unit += 0.002 * (abs(x) - abs(l)) * unit / 3
    return lev


def count_syncs(dibits):
    sh, n = 0, 0
    for d in dibits:
        sh = ((sh << 2) | int(d)) & 0xFFFFFFFFFFFF
        n += bin(sh ^ SYNC).count("1") <= 4
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wav")
    ap.add_argument("out")
    ap.add_argument("--disc", action="store_true")
    ap.add_argument("--invert", action="store_true")
    a = ap.parse_args()
    fs, v = read_wav(a.wav)
    if a.disc:
        d = v[:, 0] - np.mean(v[:, 0])
        print("discriminator WAV, %d Hz, %.1f s" % (fs, len(d) / fs))
    else:
        iq = v[:, 0] + 1j * v[:, 1]
        d, offset, active = discriminate(iq, fs)
        print("I/Q WAV, %d Hz, %.1f s; carrier offset %+.0f Hz; signal present %.0f %%" %
              (fs, len(iq) / fs, offset, 100 * active))
    s = symbols(d, fs)
    if a.invert:
        s = -s
    lev = slice_levels(s)
    dibits = np.array([LEVEL_TO_DIBIT[int(l)] for l in lev], dtype=np.uint8)
    inv = np.array([LEVEL_TO_DIBIT[-int(l)] for l in lev], dtype=np.uint8)
    n, ni = count_syncs(dibits), count_syncs(inv)
    print("%d symbols; frame syncs: %d, %d with the polarity flipped" % (len(dibits), n, ni))
    if ni > n:
        print("using the flipped polarity")
        dibits = inv
    dibits.tofile(a.out)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
