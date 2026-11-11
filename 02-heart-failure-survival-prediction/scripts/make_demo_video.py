"""Render docs/demo.mp4: a captioned, narration-free walkthrough of this project.

It RUNS the commands listed in SPEC (tests, pipeline, evaluator problem), captures the real
output, and renders a terminal-style video plus the generated figures. Re-run it after you
have run the pipeline on real data to get a video with your real results:

    python scripts/make_demo_video.py

Requires: pillow, and ffmpeg on PATH. Output: docs/demo.mp4 and docs/video_thumbnail.png.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SPEC = {'scenes': [{'type': 'card', 'title': 'Heart Failure Survival and Risk Prediction', 'lines': ['Kaplan-Meier, Cox regression, leakage-safe prediction', 'Python: lifelines, scikit-learn, pandas', 'Auto-generated walkthrough (no narration)'], 'seconds': 5, 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.'}, {'type': 'terminal', 'cmd': 'python -m pytest -q', 'caption': 'Step 1: run the test-suite', 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.', 'max_lines': 14}, {'type': 'terminal', 'cmd': 'python scripts/run_pipeline.py --synthetic', 'caption': 'Step 2: run the pipeline (synthetic cohort, no download)', 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.', 'max_lines': 22, 'hold': 4}, {'type': 'image', 'path': 'results/figures/km_ef.png', 'caption': 'Kaplan-Meier by ejection fraction (synthetic)', 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.', 'seconds': 5}, {'type': 'image', 'path': 'results/figures/calibration.png', 'caption': 'Calibration of the classifier (synthetic)', 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.', 'seconds': 5}, {'type': 'terminal', 'cmd': 'python problems/problem_01_solution.py', 'caption': 'Step 3: the evaluator problem (a mistake an AI model would likely make)', 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.', 'max_lines': 16, 'hold': 4}, {'type': 'card', 'title': 'Limitations and next steps', 'lines': ['Real run: python scripts/run_pipeline.py --csv data/heart_failure_clinical_records_dataset.csv', '299 patients, one cohort, no external validation', 'Re-render this video after your real run: python scripts/make_demo_video.py'], 'seconds': 6, 'banner': 'SYNTHETIC DEMO DATA: numbers in this video are NOT real findings. Run the pipeline on the real dataset for results.'}]}

W, H, FPS = 1280, 720, 10
BG, FG, DIM, ACC, WARN = (18, 18, 24), (225, 225, 230), (140, 140, 155), (110, 200, 255), (255, 190, 90)
MONO = ["/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", "/System/Library/Fonts/Menlo.ttc",
        "C:/Windows/Fonts/consola.ttf", "DejaVuSansMono.ttf"]
SANS = ["/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/System/Library/Fonts/Helvetica.ttc",
        "C:/Windows/Fonts/arial.ttf", "DejaVuSans.ttf"]
SANSB = ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", "C:/Windows/Fonts/arialbd.ttf",
         "DejaVuSans-Bold.ttf"]


def font(cands, size):
    for c in cands:
        try:
            return ImageFont.truetype(c, size)
        except OSError:
            continue
    return ImageFont.load_default()


F_MONO, F_SANS, F_BOLD = font(MONO, 17), font(SANS, 22), font(SANSB, 38)
F_SMALL = font(SANS, 18)
CW = F_MONO.getlength("M")
LINE_H = 23
COLS = int((W - 70) / CW)
TOP, BOTTOM = 70, 80
ROWS = (H - TOP - BOTTOM - 20) // LINE_H


class Writer:
    def __init__(self, tmp: Path):
        self.tmp, self.n = tmp, 0

    def add(self, img: Image.Image, count: int = 1):
        for _ in range(count):
            img.save(self.tmp / f"f_{self.n:05d}.png"); self.n += 1


def base(banner=None, caption=None):
    im = Image.new("RGB", (W, H), BG); d = ImageDraw.Draw(im)
    if banner:
        d.rectangle([0, 0, W, 42], fill=(60, 45, 15)); d.text((20, 9), banner, font=F_SMALL, fill=WARN)
    if caption:
        d.rectangle([0, H - BOTTOM, W, H], fill=(28, 28, 38)); d.text((24, H - BOTTOM + 24), caption, font=F_SANS, fill=FG)
    return im, d


def wrap(text, width):
    out = []
    for ln in text.splitlines() or [""]:
        ln = ln.expandtabs(4).rstrip()
        out.append(ln if len(ln) <= width else ln[: width - 1] + "…")
    return out


def card(w: Writer, title, lines, seconds=4, banner=None):
    im, d = base(banner)
    d.text((70, 150), title, font=F_BOLD, fill=ACC)
    y = 250
    for ln in lines:
        d.text((70, y), ln, font=F_SANS, fill=FG); y += 38
    w.add(im, int(seconds * FPS))
    return im


def terminal(w: Writer, cmd, caption, banner=None, max_lines=22, hold=3.0, timeout=1800):
    run = cmd.replace("python ", f'"{sys.executable}" ', 1) if cmd.startswith("python ") else cmd
    p = subprocess.run(run, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    lines = wrap((p.stdout + p.stderr).rstrip(), COLS)
    if len(lines) > max_lines:
        lines = [f"[... {len(lines) - max_lines + 1} earlier lines omitted ...]"] + lines[-(max_lines - 1):]
    prompt = "$ "
    def frame(typed, shown):
        im, d = base(banner, caption)
        d.text((35, TOP), prompt + typed + ("_" if shown == 0 else ""), font=F_MONO, fill=ACC)
        vis = lines[:shown][-ROWS:]
        for i, ln in enumerate(vis):
            d.text((35, TOP + (i + 1) * LINE_H), ln, font=F_MONO, fill=FG)
        return im
    for k in range(0, len(cmd) + 1, 3):
        w.add(frame(cmd[:k], 0))
    w.add(frame(cmd, 0), 3)
    for s in range(0, len(lines) + 1, 4):
        w.add(frame(cmd, s))
    w.add(frame(cmd, len(lines)), int(hold * FPS))
    if p.returncode != 0:
        print(f"WARNING: '{cmd}' exited with code {p.returncode}", file=sys.stderr)


def image(w: Writer, path, caption, banner=None, seconds=5):
    pth = ROOT / path
    if not pth.exists():
        print(f"WARNING: missing {path}, scene skipped", file=sys.stderr); return
    im, d = base(banner, caption)
    fig = Image.open(pth).convert("RGB")
    box = (W - 40, H - BOTTOM - (50 if banner else 10) - 20)
    fig.thumbnail(box)
    bg = Image.new("RGB", fig.size, (255, 255, 255)); bg.paste(fig)
    im.paste(bg, ((W - fig.width) // 2, (48 if banner else 10) + (box[1] - fig.height) // 2))
    w.add(im, int(seconds * FPS))


def main():
    if not shutil.which("ffmpeg"):
        sys.exit("ffmpeg not found on PATH; install it (e.g. apt install ffmpeg / brew install ffmpeg)")
    docs = ROOT / "docs"; docs.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        w = Writer(Path(t)); thumb = None
        for sc in SPEC["scenes"]:
            kind = sc["type"]
            if kind == "card":
                im = card(w, sc["title"], sc["lines"], sc.get("seconds", 4), sc.get("banner"))
                thumb = thumb or im
            elif kind == "terminal":
                terminal(w, sc["cmd"], sc["caption"], sc.get("banner"), sc.get("max_lines", 22), sc.get("hold", 3.0))
            elif kind == "image":
                image(w, sc["path"], sc["caption"], sc.get("banner"), sc.get("seconds", 5))
        thumb.save(docs / "video_thumbnail.png")
        out = docs / "demo.mp4"
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", f"{t}/f_%05d.png",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "26", "-movflags", "+faststart", str(out)],
                       check=True)
        print(f"wrote {out} ({w.n / FPS:.0f} s, {out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
