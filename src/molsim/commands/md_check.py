#import warnings # MDAnalysis throws deprecationwarning on xdrlib which is slated for removal in Python 3.13
#warnings.filterwarnings("ignore", message=".*'xdrlib' is deprecated.*", category=DeprecationWarning)

import argparse
import MDAnalysis as mda
from collections import Counter
from molsim.utils.files import getFileNames
import os

#import xdrlib
#import sys
#import numpy as np
#import re
#from math import ceil

def register(subparsers):
    parser = subparsers.add_parser('md_check',
                                   help='Print information on a `.gro` or `.pdb` file.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='md_check.py',
                                     description='Print information on a `.gro` or `.pdb` file.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('coordfile',
                        help='the coordinate file')
    return parser


def validateArguments(args):
    SUPPORTED_FILE_EXTENSIONS = [".gro", ".pdb"]
    check = False
    for ext in SUPPORTED_FILE_EXTENSIONS:
        if args.coordfile.endswith(ext):
            directory, filename = os.path.split(args.coordfile)
            check = True
            assert filename in getFileNames(ext, directory), f"Could not find {args.coordfile}"
            break
    assert check, f"`{ext}` files are not supported"
    return args


def main(args):
    args = validateArguments(args)
    """
    """
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASX", "GLX",
        "XAA", "UNK"
    ]
    Nav = 6.02214076e23  # [mol^{-1}]
    conc_water = 55.55  # [mol/L]
    ions = ["CL", "NA", "CLA", "SOD"]
    solvent = ["SOL", "H2O", "WAT", "TIP3", "TIP2"]

    print(f"Characterizing {args.coordfile}")

    ###########
    # Parsing #
    ###########
    # Create the Universe object
    u = mda.Universe(args.coordfile)
    mol_counter = Counter()
    mol_counter.update([residue.resname for residue in u.residues])
    ions_present = [ion for ion in mol_counter.keys() if ion in ions]
    solvent_pres = None
    for tag in mol_counter:
        if tag in solvent:
            solvent_pres = tag
            break

    ##############
    # Outputting #
    ##############
    print("Different kind of molecules present:")
    for val, count in mol_counter.items():
        print(f"\t{val:4s} {count}")
    dims = u.dimensions[0:3]
    cell_volume = dims[0]*dims[1]*dims[2]
    print("Box dimensions in Angstroms (xyz):")
    print("\t%-8.3f %-8.3f %-8.3f" % tuple(dims))
    # Ion concentrations
    if len(ions_present) != 0:
        assert solvent is not None, "Couldn't figure out which is the solvent"
        print("Osmolarity of ions (based on the entire box):")
        for ion in ions_present:
            print("\t%-4s %6.4f" % (ion, mol_counter[ion]/Nav/(cell_volume*1e-27)))
        if solvent_pres is not None:
            print("Molar ratio of ions to solvent:")
            for ion in ions_present:
                print("\t%-4s %11e" % (ion, float(mol_counter[ion])/float(mol_counter[solvent_pres])))
            print("Osmolarity of ions (based on molar ratio at standard conditions):")
            for ion in ions_present:
                conc = float(mol_counter[ion])/Nav/(float(mol_counter[solvent_pres])/Nav/conc_water)
                print("\t%-4s %6.4f" % (ion, conc))
    else:
        print(f"There are no recognised ions inside `{args.coordfile}`.")




if __name__ == "__main__":
    args = parseArguments()
    main(args)
