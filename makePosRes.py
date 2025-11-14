import MDAnalysis as mda
import argparse
import os
import sys

def getFileNames(extension=None, path=None, include_hidden=False):
    """
    Returns a list containing the filenames in the working directory or somewhere else
    """
    if path is None:  # Use current dir
        path = os.getcwd()
    else:
        assert os.path.exists(path)
    files_in_directory = os.listdir(path)
    # Filter out directories from the list
    if extension is None:
        filenames = [file for file in files_in_directory \
                     if os.path.isfile(os.path.join(path, file))]
    else:
        filenames = [file for file in files_in_directory \
                     if os.path.isfile(os.path.join(path, file)) \
                     and file.endswith(extension)]
    if not include_hidden:
        for index in range(len(filenames))[::-1]:
            if filenames[index].startswith('.'):
                del filenames[index]
    filenames.sort()
    return filenames


def parseArguments():
    parser = argparse.ArgumentParser(prog='makePosRes.py',
                                     description='Provide pdb or gro file for which to generate `posre.itp`. `posre.itp.` is a file that has force constants for all elements that are not `H`.',
                                     epilog='')

    parser.add_argument('filename1')           # positional argument
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

    args = parser.parse_args()
    assert args.output not in getFileNames() or args.force, f"{args.filename1} already exists, use -f flag to overwrite"

    assert '/' not in getFileNames(), "This script does not work for input not inside the working dir."
    return args

LIGHT_ELEMENTS = ['H']

args = parseArguments()

prefix_list = []
prefix_list.append("""; In this topology include file, you will find position restraint
; entries for all the heavy atoms that were inside `""")
prefix_list.append(f"{args.filename1}`.\n")
prefix_list.append("""
[ position_restraints ]
; atom  type      fx      fy      fz
""")


try:
    u = mda.Universe(args.filename1)
except Exception as e:
    sys.exit(f"Unable to parse {args.filename1}.")

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


