import argparse

def main():
    parser = argparse.ArgumentParser(prog="molsim")
    subparsers = parser.add_subparsers(dest="command", required=True) #parser.parse_args now returns a `args` with a `command` attribute

    # convert
    from commands.xvg_min import register as register_convert
    register_convert(subparsers)


    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
