#!/usr/bin/env python3
"""Supplementary Fig. 4 composite: the three Fig. 1b variants in one row, each with its own legend (AUPRCs differ) and a
label naming the variant, since the source panels share one title. Re-run after plot_sensfdr_curve.py."""
import os
from PIL import Image, ImageDraw, ImageFont

from paths import V2 as D
from paths import bold_font
PANELS = [("a", "fig1b_pr_curve_allvar.png", "All variants scored (not only credible-set members)"),
          ("b", "fig1b_pr_curve_inCS_withpal.png", "Credible-set members, A/T and C/G variants kept"),
          ("c", "fig1b_pr_curve_allvar_withpal.png", "All variants scored, A/T and C/G variants kept")]
HEAD, GAP, PAD = 170, 60, 40
ims = []
for _, f, _ in PANELS:
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i); ims.append(bg.convert("RGB"))
W = sum(i.size[0] for i in ims) + GAP * (len(ims) - 1) + 2 * PAD     # one row: easier to compare side by side
H = max(i.size[1] for i in ims) + HEAD + 2 * PAD
fig = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(fig)
fb = bold_font(120)
ft = bold_font(68)
x = PAD
for (lab, _, txt), im in zip(PANELS, ims):
    d.text((x, PAD), lab, font=fb, fill="black")
    d.text((x + 120, PAD + 38), txt, font=ft, fill="#3A3A3A")
    fig.paste(im, (x, PAD + HEAD))
    x += im.size[0] + GAP
fig.save(os.path.join(D, "supp_fig4_pr_variants.png"), dpi=(300, 300))
print("wrote supp_fig4_pr_variants.png", fig.size)
