#!/usr/bin/env python3
"""Supplementary Fig. 8 composite: kriging at the Fig. 2b locus (a EUR, b EAS; plot_supp_kriging.py) on top, and the in-CS
precision-recall of Fig. 1b with the high-s loci removed (c; plot_sensfdr_curve.py --s-flag 0.01) below, at C_FRAC of
the width so its fonts match a/b. Panel c answers a co-author's round-3 comment: the high-PIP dip in SuSiEx's curve is the two calls at this locus, which has
the largest per-ancestry s of all 250 (B70, EAS 0.076). Re-run after both plotting scripts, then assemble_supp.py."""
import os
from PIL import Image, ImageDraw, ImageFont

from paths import V2 as D
from paths import bold_font
TOP, BOT, OUT = "supp_kriging_chr12_29898169.png", "supp_pr_inCS_sflag.png", "supp_fig8_kriging_pr.png"
GAP, PAD = 60, 40
C_FRAC = 0.72   # c at full width set its fonts ~1.4x a/b's; this matches them


def rgb(f):
    i = Image.open(os.path.join(D, f)).convert("RGBA")
    bg = Image.new("RGBA", i.size, "white"); bg.alpha_composite(i)
    return bg.convert("RGB")


top, bot = rgb(TOP), rgb(BOT)
w = top.size[0]
wc = round(w * C_FRAC)
bot = bot.resize((wc, round(bot.size[1] * wc / bot.size[0])), Image.LANCZOS)
fig = Image.new("RGB", (w + 2 * PAD, top.size[1] + GAP + bot.size[1] + 2 * PAD), "white")
fig.paste(top, (PAD, PAD))
y = PAD + top.size[1] + GAP
fig.paste(bot, (PAD, y))
# same size and inset as the a/b labels drawn by plot_supp_kriging.py (16 pt at 300 dpi)
fb = bold_font(70)
ImageDraw.Draw(fig).text((PAD + 4, y - 10), "c", font=fb, fill="#1a1a1a")
fig.save(os.path.join(D, OUT), dpi=(300, 300))
print("wrote", OUT, fig.size)
