#!/usr/bin/env python3
"""Fig 2 composite: a-c in one row, d below. a-c are the narrow locus panels (plot_fig2_locus.py, 6 in wide); d is plot_fig2d.py at its
own size (figures/v2/fig2d_top13.png), top-aligned, so a-c keep most of the width. Shorter panels are padded with
white at the bottom, never stretched. Swap 2c with --c (e.g. the chr1 alternate)."""
import argparse
import os
from PIL import Image, ImageDraw, ImageFont

from paths import V2 as D
from paths import bold_font
ap = argparse.ArgumentParser()
ap.add_argument("--c", default="fig2c_chr3_36834083_36998583.png")
a = ap.parse_args()
files = ["fig2a_chr12_110129434_111489434.png", "fig2b_chr12_29898169_29948269.png", a.c, "fig2d_top13_wide.png"]   # 2d transposed (plot_fig2d.py --transpose), one row under a-c
ims = []
for f in files:
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i); ims.append(bg.convert("RGB"))
GAP, TOP, SIDE = 60, 150, 40
# two rows: a-c across the top at full width, d below on its own, left-aligned at its native size
top, d4 = ims[:3], ims[3]
H1 = max(i.size[1] for i in top)
W = SIDE * 2 + sum(i.size[0] for i in top) + GAP * 2
fig = Image.new("RGB", (W, TOP + H1 + GAP + TOP + d4.size[1] + SIDE), "white")
d = ImageDraw.Draw(fig)
font = bold_font(130)
x = SIDE
for lab, im in zip("abc", top):
    fig.paste(im, (x, TOP))
    d.text((x, 5), lab, font=font, fill="black")
    x += im.size[0] + GAP
y2 = TOP + H1 + GAP
d.text((SIDE, y2 + 5), "d", font=font, fill="black")
fig.paste(d4, (SIDE, y2 + TOP))
fig.save(os.path.join(D, "fig2.png"), dpi=(300, 300))
print("wrote fig2.png", fig.size, "2c =", a.c)
