import MDAnalysis as mda
import argparse
import sys


from molsim.utils import getFileNames


def register(subparsers):
    parser = subparsers.add_parser('make_posres',
                                   help='Provide pdb or gro file for which to generate `posre.itp`. `posre.itp.` is a file that has force constants for all elements that are not `H`.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='make_posres.py',
                                     description='Provide pdb or gro file for which to generate `posre.itp`. `posre.itp.` is a file that has force constants for all elements that are not `H`.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('coordfile',
                        help='the coordinate file')           # positional argument
    parser.add_argument("-o",
                        "--output",
                        default="posre.itp",
                        type=str,
                        nargs='?',
                        help='name of the output which defaults to `posre.itp`')
    parser.add_argument("-f",
                        "--force",
                        action='store_true',
                        help="set to overwrite existing files")
    return parser


def validateArguments(args):
    assert args.output not in getFileNames() or args.force, f"{args.coordfile} already exists, use -f flag to overwrite"

    return args


def main(args):
    args = validateArguments(args)


    LIGHT_ELEMENTS = ['H']


    prefix_list = []
    prefix_list.append("""; In this topology include file, you will find position restraint
; entries for all the heavy atoms that were inside `""")
    prefix_list.append(f"{args.coordfile}`.\n")
    prefix_list.append("""
[ position_restraints ]
; atom  type    fx    fy    fz
""")


    try:
        u = mda.Universe(args.coordfile)
    except Exception:
        sys.exit(f"Unable to parse {args.coordfile}.")

    with open(args.output, 'w') as file:
        for prefix in prefix_list:
            file.write(prefix)

        for atom in u.atoms:
            Heavy = True
            if hasattr(atom, 'element') and atom.element not in LIGHT_ELEMENTS:
                Heavy = False
            else:
                for name in LIGHT_ELEMENTS:
                    if atom.name.startswith(name):
                        Heavy = False
                if Heavy:
                    file.write(f"{atom.id:6d}{1:6d}{1000:6}{1000:6}{1000:6}\n")


if __name__ == "__main__":
    args = parseArguments()
    main(args)
