"""Our IMBE encoder against the reference encoder on real speech (the THCHS-30 sentences of dmr-vocoder): each
sentence encoded by ours (bin/imbe_enc_cli.exe) and by the reference (bin/imbe_ref.exe enc, c/build.sh ref), both
decoded by the reference decoder and by ours; PESQ against the original, the output level, and the parameters side by
side (bin/imbe_dump.exe: pitch agreement, voiced harmonics, the mean log2 amplitude of voiced and unvoiced ones).

    py py/enc_test.py [sentence ...]
"""

import glob
import os
import subprocess
import sys

import numpy as np
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
BIN = os.path.join(ROOT, "bin")
SPEECH = os.path.join(ROOT, "..", "vocoder", "data", "thchs30")
PESQ = os.path.join(ROOT, "..", "vocoder", "bin", "pesq_nb.exe")
OUT = os.path.join(ROOT, "out", "enc")


def run(*args):
    r = subprocess.run(list(args), capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("%s: %s" % (args[0], r.stderr or r.stdout))
    return r.stdout


def wav_samples(path):
    return np.fromfile(path, dtype="<i2", offset=44)


def pesq(ref_raw, deg_raw):
    return float(run(PESQ, ref_raw, deg_raw).split()[-1])


def dump(frames):
    return np.array([[float(v) for v in l.split()] for l in run(os.path.join(BIN, "imbe_dump.exe"), frames).splitlines()])


def main():
    os.makedirs(OUT, exist_ok=True)
    names = sys.argv[1:] or [os.path.splitext(os.path.basename(f))[0] for f in sorted(glob.glob(os.path.join(SPEECH, "*.wav")))]
    tot = {}
    for name in names:
        x = wav_samples(os.path.join(SPEECH, name + ".wav")).astype(np.float64)
        x8 = np.clip(resample_poly(x, 1, 2), -32768, 32767).astype("<i2")
        orig = os.path.join(OUT, name + ".raw")
        x8.tofile(orig)
        enc = {}
        for who, exe, args in (("ours", "imbe_enc_cli.exe", []), ("ref", "imbe_ref.exe", ["enc"])):
            enc[who] = os.path.join(OUT, "%s_%s.imbe" % (name, who))
            run(os.path.join(BIN, exe), *args, orig, enc[who])
        line = "%-8s" % name
        for e in ("ours", "ref"):
            for d, exe, args in (("ours", "imbe_cli.exe", []), ("ref", "imbe_ref.exe", ["dec"])):
                wav = os.path.join(OUT, "%s_%s_%s.wav" % (name, e, d))
                run(os.path.join(BIN, exe), *args, enc[e], wav)
                y = wav_samples(wav)
                raw = wav[:-4] + ".raw"
                y.tofile(raw)
                p = pesq(orig, raw)
                lv = 20 * np.log10(np.sqrt(np.mean(y.astype(float) ** 2)) / np.sqrt(np.mean(x8.astype(float) ** 2)))
                tot.setdefault((e, d), []).append(p)
                line += "  enc %s/dec %s %.3f %+5.1f dB" % (e, d, p, lv)
        a, b = dump(enc["ours"]), dump(enc["ref"])
        n = min(len(a), len(b))
        best = None
        for sh in range(-3, 4):  # the encoders' delays may differ
            i0, j0 = max(0, sh), max(0, -sh)
            m = n - abs(sh)
            agree = np.mean(np.abs(a[i0:i0 + m, 0] - b[j0:j0 + m, 0]) <= 4)
            if best is None or agree > best[0]:
                best = (agree, sh)
        agree, sh = best
        i0, j0 = max(0, sh), max(0, -sh)
        m = n - abs(sh)
        A, B = a[i0:i0 + m], b[j0:j0 + m]
        both_v = (A[:, 2] > 0) & (B[:, 2] > 0)
        both_u = (A[:, 2] < A[:, 1]) & (B[:, 2] < B[:, 1])
        print(line)
        print("          pitch within 4 steps %.0f%% (shift %d); voiced harmonics ours %.1f ref %.1f; "
              "mean log2 voiced ours-ref %+.2f, unvoiced %+.2f" %
              (100 * agree, sh, A[:, 2].mean(), B[:, 2].mean(), (A[both_v, 3] - B[both_v, 3]).mean(),
               (A[both_u, 4] - B[both_u, 4]).mean()))
    print("mean " + "  ".join("enc %s/dec %s %.3f" % (e, d, np.mean(v)) for (e, d), v in tot.items()))


if __name__ == "__main__":
    main()
