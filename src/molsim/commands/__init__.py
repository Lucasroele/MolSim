# Register the available subcommands
from .xvg_min import register as register_xvg_min
from .md_check import register as register_md_check
from .make_posres import register as register_make_posres
from .plot_xvg import register as register_plot_xvg

__all__ = ['register_xvg_min', 'register_md_check', 'register_make_posres', 'register_plot_xvg']

