"""Where the v2 figures and tables are read from and written to.

Defaults to figures/v2/ next to this file. Set V2_OUT to rebuild everything into another
directory instead (regenerate_v2.sh does this to check a rebuild against the committed files).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.environ.get("V2_OUT", os.path.join(HERE, "figures", "v2")))


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
