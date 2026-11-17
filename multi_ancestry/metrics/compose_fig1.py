#!/usr/bin/env python3
"""Fig 1 composite (Qin): a = study-design schematic (fig1A.png, Qin's, flattened onto white) on the left,
b = in-CS PR curve top right, c = exact lead recovery bottom right. No resampling up: the right column is
set to 1c's native width and 1b is scaled down to it; the left column is scaled to the right column's height."""
import os
from PIL import Image, ImageDraw, ImageFont

from paths import V2 as D
GAP, PAD, LPAD = 70, 40, 130   # LPAD: room left of panel a for its letter


def white(f):
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i)
    return bg.convert("RGB")


a, b, c = white("fig1A.png"), white("fig1b_pr_curve_inCS.png"), white("fig1c_recovery.png")
Wr = c.size[0]
b = b.resize((Wr, round(b.size[1] * Wr / b.size[0])), Image.LANCZOS)
H = b.size[1] + GAP + c.size[1]
a = a.resize((round(a.size[0] * H / a.size[1]), H), Image.LANCZOS)
W = LPAD + a.size[0] + GAP + Wr + PAD
fig = Image.new("RGB", (W, H + 2 * PAD), "white")
xa, xr = LPAD, LPAD + a.size[0] + GAP
fig.paste(a, (xa, PAD)); fig.paste(b, (xr, PAD)); fig.paste(c, (xr, PAD + b.size[1] + GAP))
d = ImageDraw.Draw(fig)
f = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 110)
for lab, (x, y) in {"a": (10, 0), "b": (xr - 60, 0), "c": (xr - 60, PAD + b.size[1] + GAP - 20)}.items():
    d.text((max(x, 0), y), lab, font=f, fill="black")
fig.save(os.path.join(D, "fig1.png"), dpi=(300, 300))
print("wrote fig1.png", fig.size)
