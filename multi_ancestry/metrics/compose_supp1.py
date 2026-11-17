#!/usr/bin/env python3
"""Supplementary Fig. 1 composite: a = the hand-drawn data flowchart ("supp fig1A.png", hand-made, flattened onto white) at its
native size, which also sets the canvas; b = the variants each method fitted (plot_variant_domains.py --v2 --arms-only,
supp_fig1b_arm_domains.png), scaled (LANCZOS) to 2,160 px wide and pasted into the flowchart's empty top-right corner,
20 px from the top and right edges. Panel letters Arial Bold 96. Re-run after plot_variant_domains.py, then
assemble_supp.py. Panel b's empty top band is trimmed before pasting (2026-10-11).
  python3 compose_supp1.py [--dir figures/v2]
"""
import argparse
from paths import V2
from paths import bold_font
import os
from PIL import Image, ImageDraw, ImageFont, ImageOps

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", default=V2,
                help="reads both panels from and writes supp_fig1_design_domains.png to this directory")
D = ap.parse_args().dir
MARGIN = 60
WB, EDGE = 2080, 20   # panel b width (2160 until 2026-10-11, when its legend crowded the flowchart); margin from top and right edges


def white(f):
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i)
    return bg.convert("RGB")


a, b = white("supp fig1A.png"), white("supp_fig1b_arm_domains.png")
# drop panel b's empty top band (2026-10-11: b's legend crowded the flowchart below), keeping a 10 px margin
top = Image.eval(b.convert("L"), lambda v: 255 if v < 250 else 0).getbbox()[1]
b = b.crop((0, max(0, top - 10), b.size[0], b.size[1]))
b = b.resize((WB, int(b.size[1] * WB / b.size[0])), Image.LANCZOS)
fig = a.copy()
xb = fig.size[0] - EDGE - WB
fig.paste(b, (xb, EDGE))
d = ImageDraw.Draw(fig)
f = bold_font(96)
for lab, (x, y) in {"a": (12, 0), "b": (xb - 70, 0)}.items():
    d.text((x, y), lab, font=f, fill="black")
# white margin all round, so the flowchart's frames do not touch the figure edge
fig = ImageOps.expand(fig, border=MARGIN, fill="white")
fig.save(os.path.join(D, "supp_fig1_design_domains.png"), dpi=(300, 300))
print("wrote supp_fig1_design_domains.png", fig.size)
