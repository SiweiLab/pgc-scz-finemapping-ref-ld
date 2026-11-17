"""Where the v2 figures and tables are read from and written to.

Defaults to figures/v2/ next to this file. Set V2_OUT to rebuild everything into another
directory instead (regenerate_v2.sh does this to check a rebuild against the committed files).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.environ.get("V2_OUT", os.path.join(HERE, "figures", "v2")))
