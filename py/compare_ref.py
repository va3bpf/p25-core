"""Compare our IMBE decoder's output with a reference decoder's on the same frames (both WAV, 8 kHz):
level, delay, log-spectral distance, and PESQ of ours with the reference as the reference signal
(dmr-vocoder's bin/pesq_nb.exe, next to this repository in OpenDMR).

    py py/compare_ref.py ours.wav ref.wav
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PESQ = os.path.join(HERE, "..", "..", "vocoder", "bin", "pesq_nb.exe")


def read(path):
    return np.fromfile(path, dtype="<i2", offset=44).astype(np.float64)


def db(x):
    return 20 * np.log10(np.sqrt(np.mean(x ** 2)) + 1e-9)


def pesq(ref, deg):
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, "ref.raw"), os.path.join(d, "deg.raw")
        np.clip(ref, -32768, 32767).astype("<i2").tofile(a)
        np.clip(deg, -32768, 32767).astype("<i2").tofile(b)
        r = subprocess.run([PESQ, a, b], capture_output=True, text=True)
        return r.stdout.strip() or r.stderr.strip()


def main():
    ours, ref = read(sys.argv[1]), read(sys.argv[2])
    n = min(len(ours), len(ref))
    ours, ref = ours[:n], ref[:n]
    active = np.abs(ref) > 0.05 * np.max(np.abs(ref))
    print("level: ours %.1f dBFS, reference %.1f dBFS (ours %+.1f dB)" %
          (db(ours) - 90.3, db(ref) - 90.3, db(ours) - db(ref)))
    # delay: from the 5 ms energy envelopes (the waveforms differ: each decoder makes its own harmonic phases)
    env = lambda x: np.array([np.sqrt(np.mean(x[i:i + 40] ** 2)) for i in range(0, n - 40, 40)])
    eo, er = env(ours), env(ref)
    best = max(range(-12, 13), key=lambda l: np.corrcoef(np.roll(eo, -l), er)[0, 1])
    lag = best * 40
    print("delay of ours: %+.0f ms (envelope correlation %.3f)" %
          (lag / 8, np.corrcoef(np.roll(eo, -best), er)[0, 1]))
    # the spectral envelopes (smoothed over 8 bins, 200..3400 Hz) on 32 ms frames where the reference speaks,
    # the levels matched
    g = np.sqrt(np.mean(ref ** 2) / np.mean(ours ** 2))
    o = np.roll(ours * g, -lag)
    lsd = []
    sm = lambda X: np.convolve(X, np.ones(8) / 8, mode="same")
    band = slice(6, 109)
    for i in range(0, n - 256, 160):
        if np.mean(active[i:i + 256]) < 0.5:
            continue
        w = np.hanning(256)
        A = sm(10 * np.log10(np.abs(np.fft.rfft(ref[i:i + 256] * w)) ** 2 + 1e-3))
        B = sm(10 * np.log10(np.abs(np.fft.rfft(o[i:i + 256] * w)) ** 2 + 1e-3))
        lsd.append(np.sqrt(np.mean((A[band] - B[band]) ** 2)))
    print("spectral envelope distance (speech frames): mean %.1f dB over %d frames" % (np.mean(lsd), len(lsd)))
    if os.path.exists(PESQ):
        print("PESQ of ours against the reference:", pesq(ref, o))
        print("PESQ of the reference against itself:", pesq(ref, ref))
    else:
        print("(no PESQ: %s not built)" % PESQ)


if __name__ == "__main__":
    main()
