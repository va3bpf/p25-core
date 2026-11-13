"""Our IMBE decoder against a reference decoder on real speech: the THCHS-30 sentences of dmr-vocoder
(vocoder/data/thchs30, next to this repository in OpenDMR) encoded with the reference encoder (bin/imbe_ref.exe,
c/build.sh ref), the frames decoded by ours (bin/imbe_cli.exe) and by the reference, PESQ of each against the
original speech.

    py py/speech_test.py
"""

import glob
import os
import subprocess

import numpy as np
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
SPEECH = os.path.join(ROOT, "..", "vocoder", "data", "thchs30")
PESQ = os.path.join(ROOT, "..", "vocoder", "bin", "pesq_nb.exe")
OUT = os.path.join(ROOT, "out", "speech")


def run(*args):
    r = subprocess.run(list(args), capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("%s: %s" % (args[0], r.stderr or r.stdout))
    return r.stdout


def wav_samples(path):
    return np.fromfile(path, dtype="<i2", offset=44)


def pesq(ref_raw, deg_raw):
    return float(run(PESQ, ref_raw, deg_raw).split()[-1])  # MOS-LQO (P.862.1)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(SPEECH, "*.wav"))):
        name = os.path.splitext(os.path.basename(f))[0]
        x = wav_samples(f).astype(np.float64)
        x8 = np.clip(resample_poly(x, 1, 2), -32768, 32767).astype("<i2")
        orig = os.path.join(OUT, name + ".raw")
        x8.tofile(orig)
        frames = os.path.join(OUT, name + ".imbe")
        run(os.path.join(ROOT, "bin", "imbe_ref.exe"), "enc", orig, frames)
        res = {}
        for who, exe, args in (("ours", "imbe_cli.exe", []), ("ref", "imbe_ref.exe", ["dec"])):
            wav = os.path.join(OUT, "%s_%s.wav" % (name, who))
            run(os.path.join(ROOT, "bin", exe), *args, frames, wav)
            raw = os.path.join(OUT, "%s_%s.raw" % (name, who))
            y = wav_samples(wav)
            y.tofile(raw)
            res[who] = (pesq(orig, raw), 20 * np.log10(np.sqrt(np.mean(y.astype(float) ** 2)) /
                                                         np.sqrt(np.mean(x8.astype(float) ** 2))))
        rows.append((name, res))
        print("%-9s PESQ ours %.3f  reference %.3f   level vs original: ours %+.1f dB, reference %+.1f dB" %
              (name, res["ours"][0], res["ref"][0], res["ours"][1], res["ref"][1]))
    po = np.mean([r["ours"][0] for _, r in rows])
    pr = np.mean([r["ref"][0] for _, r in rows])
    print("mean      PESQ ours %.3f  reference %.3f  (ours %+.3f)" % (po, pr, po - pr))


if __name__ == "__main__":
    main()
