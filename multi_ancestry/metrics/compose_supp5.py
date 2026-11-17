#!/usr/bin/env python3
"""Supplementary Fig. 5 composite: the L sweep under credible-set scoring (a, as Fig. 1b) above the same sweep with all
variants scored (b). Supplementary Note 4 rests on the contrast between the two. Re-run after plot_Lsweep_curve.py, then
assemble_supp.py."""
import os
from PIL import Image, ImageDraw, ImageFont

from paths import V2 as D
PANELS = [("a", "supp_Lsweep_pr_inCS.png", "Credible-set members scored (as Fig. 1b)"),
          ("b", "supp_Lsweep_pr_allvar.png", "All variants scored")]
HEAD, GAP, PAD = 190, 90, 50
ims = []
for _, f, _ in PANELS:
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i); ims.append(bg.convert("RGB"))
W = max(i.size[0] for i in ims) + 2 * PAD
H = sum(i.size[1] + HEAD for i in ims) + GAP * (len(ims) - 1) + 2 * PAD
fig = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(fig)
fb = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 140)
ft = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 80)
y = PAD
for (lab, _, txt), im in zip(PANELS, ims):
    d.text((PAD, y), lab, font=fb, fill="black")
    d.text((PAD + 140, y + 45), txt, font=ft, fill="#3A3A3A")
    fig.paste(im, (PAD, y + HEAD))
    y += HEAD + im.size[1] + GAP
fig.save(os.path.join(D, "supp_fig5_Lsweep.png"), dpi=(300, 300))
print("wrote supp_fig5_Lsweep.png", fig.size)
