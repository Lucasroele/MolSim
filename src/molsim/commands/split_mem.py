import argparse
import MDAnalysis as mda
from collections import Counter
from pathlib import Path
from molsim.utils.mda import isPhospholipid


def register(subparsers):
    parser = subparsers.add_parser('split_mem',
                                   help='Append or write groups to the .ndx that split the phospholipids based on whether they reside above or below the center of geometry of the membrane along the Z-axis.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='split_mem.py',
                                     description='Append or write groups to the .ndx that split the phospholipids based on whether they reside above or below the center of geometry of the membrane along the Z-axis.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('filename1',
                        help='the coordinate file that contains a lipid bilayer.')           # positional argument
    parser.add_argument('filename2',
                        help='the .tpr file')
    parser.add_argument('-o',
                        '--output',
                        default="index.ndx",
                        help='the ndx file that will see the Z+ and Z- groups appended.')           # positional argument
    parser.add_argument('-u',
                        '--unequal',
                        action='store_true',
                        help='set when leaflet phospholipid distribution is not identical.')
    parser.add_argument('-a',
                        '--append',
                        action="store_true",
                        help='creates a new .ndx and overwrites any .ndx present')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='set to overwrite existing files')
    parser.add_argument('-q',
                        '--quiet',
                        action='store_true',
                        help='set to suppress warnings and informational messages')
    return parser


def validateArguments(args):
    filepath1 = Path(args.filename1)
    args.filename1 = str(filepath1)
    filepath2 = Path(args.filename2)
    args.filename2 = str(filepath2)
    outpath = Path(args.output)
    if not outpath.suffix == ".ndx":
        args.output += ".ndx"
        outpath = Path(args.output)
    assert filepath1.is_file(), "Error: The file `" + filepath1.name + "` does not exist."
    assert filepath2.is_file(), "Error: The file `" + filepath2.name + "` does not exist."
    if outpath.is_file():
        assert args.force or args.append, "Error: The file `" + outpath.name + "` already exists, set `-f` to overwrite or `-a` to append."
    else:
        assert not args.append, "Error: The file `" + outpath.name + "` does not exist and you set `--append`."


    return args

def main(args):
    """
    This script creates gmx workable groups and stores them in an ndx file.
    The groups are created based on whether or not they are phospholipids and reside above or below the center of geometry of all phospholipids.

    INPUT
        .tpr and .gro
    OUTPUT
        .ndx

    python3 memSplitNDX.py in.tpr in.gro -o index.ndx
    """
    args = validateArguments(args)
    # Create the Universe object
    #u = mda.Universe(args.filename1, args.filename2)
    u = mda.Universe(args.filename2)
    residue_counter = Counter()
    residue_counter.update([residue.resname for residue in u.residues])
    resnames = [key for key, value in residue_counter.items()]
    z_plus_minus = ["Z+", "Z-"]

    # Find the phospholipids
    resids = {}  # {"<resname>": <resid>}   ###  first resid from each resname
    phospholipids = [] # ["<resname>"]   ###  resnames belonging to phospholipids
    for resname in resnames:
        for residue in u.residues:
            if residue.resname == resname:
                resids[resname] = residue.resid
                break
        molecule = u.select_atoms('resid ' + str(resids[resname]))
        if isPhospholipid(molecule):
            phospholipids.append(resname)
    assert len(phospholipids) > 0, f"No phospholipids in {args.filename2}."

    # Create two groups for every phospholipid
    cog_mem = u.select_atoms('resname ' + ' '.join([resname for resname in phospholipids])).center_of_geometry()
    phoslip_resids = {}  # {"<resname>": {"Z+": [<resid>], "Z-": [<resid>]}}
    above = Counter()
    below = Counter()
    for resname in phospholipids:
        phoslip_resids[resname] = {}
        for z in z_plus_minus:
            phoslip_resids[resname][z] = []
    for residue in u.residues:
        if residue.resname in phospholipids:
            if u.select_atoms('resid ' + str(residue.resid)).center_of_geometry()[2] > cog_mem[2]:
                phoslip_resids[residue.resname][z_plus_minus[0]].append(residue.resid)
                above.update([residue.resname])
            else:
                phoslip_resids[residue.resname][z_plus_minus[1]].append(residue.resid)
                below.update([residue.resname])

    # Make sure same number of phospholipids in both leaflets
    if not args.unequal:
        for resname in phospholipids:
            assert above[resname] == below[resname], f"Number of phospholipids in each leaflet is not equal for {resname}."

    universe_groups = {} # {"<resname>": [<atomGroup>, <atomGroup>]}
    for resname in phospholipids:
        universe_groups[resname] = {}
        for z in z_plus_minus:
            str_resids = ' '.join([str(num) for num in phoslip_resids[resname][z]])
            universe_groups[resname][z] = u.select_atoms('resid ' + str_resids)

    # Output
    if args.append:
        writemode = 'a'
    else:
        writemode = 'w'
    with mda.selections.gromacs.SelectionWriter(args.output, mode=writemode) as ndx:
        for resname in universe_groups:
            for z in z_plus_minus:
                ndx.write(universe_groups[resname][z], name=f"{resname}-" + z)
                ndx._outfile.write("\n")

    # Log
    if not args.quiet:
        print("\nCreated the following atom groups:")
        for resname in universe_groups:
            for groupname in universe_groups[resname]:
                print(f"    {resname}_{groupname}")
        if args.append:
            print("\nThey were appended to:")
        else:
            print("\nThey were written to:")
        print(f"    {args.output}")


if __name__ == '__main__':
    args = parseArguments()
    main(args)

