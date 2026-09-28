import argparse
import sys
from molsim.utils.errors import ArgumentValidationError

def main():
    parser = argparse.ArgumentParser(prog="molsim")
    subparsers = parser.add_subparsers(dest="command", required=True) #parser.parse_args now returns a `args` with a `command` attribute

    # convert
    from .commands import register_xvg_min
    register_xvg_min(subparsers)

    from .commands import register_md_check
    register_md_check(subparsers)

    from .commands import register_make_posres
    register_make_posres(subparsers)

    from .commands import register_plot_xvg
    register_plot_xvg(subparsers)

    from .commands import register_frame_times
    register_frame_times(subparsers)

    from .commands import register_top_cyclr
    register_top_cyclr(subparsers)

    from .commands import register_gen_cpep
    register_gen_cpep(subparsers)

    from .commands import register_count_mols
    register_count_mols(subparsers)

    from .commands import register_scale_coords
    register_scale_coords(subparsers)

#    from .commands import register_split_mem
#    register_split_mem(subparsers)

    from .commands import register_lipid_ndx
    register_lipid_ndx(subparsers)

    from .dev.merge_pdbs import register
    register(subparsers)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()  # 'command' is set to subcommand (e.g., xvg_min), 'func' points to the handler
    try:
        args.func(args) # Call the function associated with the chosen subcommand (main from xvg_min.py)
    except ArgumentValidationError as error:
        RED = "\033[31m"
        RESET = "\033[0m"
        print(f"molsim {args.command}: {RED}input error{RESET}:\n {error}", file=sys.stderr)
        return 2

if __name__ == "__main__":
    main()
