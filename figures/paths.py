"""Where the figures and tables are written (and the composites read their panels from).

Defaults to output/ next to this file. Set FIG_OUT to write everything to another directory
instead. The two hand-made panels live in hand_made/; make_figures.sh copies them into the output directory.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.environ.get("FIG_OUT", os.path.join(HERE, "output")))
HAND_MADE = os.path.join(HERE, "hand_made")


def bold_font(size):
    """Panel-letter font for the compose_*.py composites: Arial Bold if installed (BOLD_FONT overrides),
    else matplotlib's bundled DejaVu Sans Bold, so the composites also build off macOS."""
    from PIL import ImageFont
    for f in (os.environ.get("BOLD_FONT", ""), "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
              "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf"):
        if f and os.path.exists(f):
            return ImageFont.truetype(f, size)
    import matplotlib
    return ImageFont.truetype(os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf"), size)
