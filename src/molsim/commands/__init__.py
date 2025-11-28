# Register the available subcommands
from .xvg_min import register as register_xvg_min
from .md_check import register as register_md_check
from .make_posres import register as register_make_posres

__all__ = ['register_xvg_min', 'register_md_check', 'register_make_posres']

