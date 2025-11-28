import argparse

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

    args = parser.parse_args()  # 'command' is set to subcommand (e.g., xvg_min), 'func' points to the handler
    args.func(args)              # Call the function associated with the chosen subcommand (main from xvg_min.py)

if __name__ == "__main__":
    main()
