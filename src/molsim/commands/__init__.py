# Register the available subcommands
from .xvg_min import register as register_xvg_min
from .md_check import register as register_md_check
from .make_posres import register as register_make_posres
from .plot_xvg import register as register_plot_xvg
from .frame_times import register as register_frame_times
from .top_cyclr import register as register_top_cyclr

__all__ = ['register_xvg_min', 'register_md_check', 'register_make_posres', 'register_plot_xvg', 'register_frame_times', 'register_top_cyclr']

