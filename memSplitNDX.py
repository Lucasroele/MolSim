import warnings # MDAnalysis throws deprecationwarning on xdrlib which is slated for removal in Python 3.13
warnings.filterwarnings("ignore", message=".*'xdrlib' is deprecated.*", category=DeprecationWarning)

import xdrlib
import os
import sys
import argparse
import numpy as np
import MDAnalysis as mda
from collections import Counter


def parseArguments():
    parser = argparse.ArgumentParser(prog='memSplitNDX.py',
                                     description='Append or write groups to the .ndx that split the phospholipids based on whether they reside above or below the center of geometry of the membrane along the Z-axis.',
                                     epilog='Written by Lucas Roeleveld')

    parser.add_argument('filename1',
                        help='the coordinate file that contains a lipid bilayer.')           # positional argument
    parser.add_argument('filename2',
                        help='the .tpr file')
    parser.add_argument('-o',
                        '--output',
                        default="index.ndx",
                        help='the ndx file that will see the Z+ and Z- groups appended.')           # positional argument
    parser.add_argument('-n',
                        '--newNDX',
                        action="store_true",
                        help='creates a new .ndx and overwrites any .ndx present')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='enable this flag when leaflet phospholipid distribution is not identical.')
    args = parser.parse_args()
    return args


def isPhospholipid(atomgroup) -> bool:

    """
    Returns:
        True  -> atomgroup is a phospholipid
        False -> atomgroup is not a phospholipid (probably)

    Might not work when:
        - atomgroup contains multiple phosphates

    """
    def hasPhosGroup(atomgroup) -> np.ndarray | None:
        """
        Returns:
            None -> atomgroup does not contain a phosphate
            list -> atomgroup contains a phosphate (list contains indices of the phosphate and its oxygens)
                [<P_index> <O_index> <O_index> <O_index> <O_index>]
        """
        Ps = atomgroup.select_atoms('element P')
        if not Ps:
            return None
        else:
            for P in Ps:
                Os = []
                for bond in P.bonds:
                    if bond.partner(P).element == "O":
                        Os.append(bond.partner(P).index)
                if len(Os) == 4:
                    has_phos = []
                    has_phos.append(P.index)
                    has_phos.extend(Os)
                    return np.array(has_phos)
        return None


    def hasGlycerolAttached(atomgroup, phos_indices) -> (None, None) | (np.ndarray, np.ndarray):
        """
        Check if atomindex in molecule is an O-C-C or O-C-2C
        Returns:
            None, None
            np.array, np.array -> indices of the glycerol, indices of the carbonyl carbons attached to it
        """
        def getEsterCarbonyl(bridgeO, alcoholC):
            """
            bridgeO  = mda.Atom
            alcoholC = mda.Atom
            Returns:
                list of length 2 with indices of the carbon and oxygen of the carbonyl attached to the alcohol
            """
            carbonyl_indices = [] # [<Cindex>, <Oindex>]
            other_atoms = 0
            for bond in bridgeO.bonds:
                if bond.partner(bridgeO) == alcoholC:
                    continue
                elif bond.partner(bridgeO).element == "C":
                    carbonyl_indices.append(bond.partner(bridgeO).index)
                    for bond2 in bond.partner(bridgeO).bonds:
                        if bond2.partner(bond.partner(bridgeO)) == bridgeO:
                            continue
                        elif bond2.partner(bond.partner(bridgeO)).element == "O":
                            carbonyl_indices.append(bond2.partner(bond.partner(bridgeO)).index)
                        elif other_atoms > 1:
                            return None
                        else:
                            other_atoms += 1
                else:
                    return None
            return carbonyl_indices

        def isGlycerolCarbon(carbontype, atom):
            """
            atom = mda.Atom (the potential glycerol carbon)
            carbontype = 2:
                atom should be terminal glycerol carbon
            carbontype = 3:
                atom should be interstitial glycerol carbon
            """
            CandO = Counter()
            for bond in atom.bonds:
                CandO.update([bond.partner(atom).element])
            if not ("C" in CandO and "O" in CandO):
                return False
            if carbontype == 2:
                if CandO["C"] == 1 and CandO["O"] == 1:
                    return True
                else:
                    return False
            elif carbontype == 3:
                if CandO["C"] == 2 and CandO["O"] == 1:
                    return True
                else:
                    return False

        glycerol_found = False
        carbonylC_indices = []
        current_atom = None # the index of the atom the decider is at
        phosmol = atomgroup.select_atoms('index ' +  ' '.join([str(num) for num in phos_indices]))
        for oxygen in phosmol:
            if glycerol_found:
                break
            glyc_indices = []
            if oxygen.element != "O" or len(oxygen.bonds) == 1:
                continue
            if oxygen.bonds[0].partner(oxygen).element == "P":
                current_atom = oxygen.bonds[1].partner(oxygen)
            else:
                current_atom = oxygen.bonds[0].partner(oxygen)
            glyc_indices.append(current_atom.index) # add first glycerol C
            if isGlycerolCarbon(2,current_atom):
                for bond in current_atom.bonds:
                    if bond.partner(current_atom).element == "C":
                        prev_atom = current_atom
                        current_atom = bond.partner(current_atom)
                        glyc_indices.append(current_atom.index) # add second glycerol C
                        break
                if isGlycerolCarbon(3,current_atom):
                    for bond in current_atom.bonds:
                        if prev_atom == bond.partner(current_atom) or bond.partner(current_atom).element == "H":
                            continue
                        elif bond.partner(current_atom).element == "C":
                            next_atom = bond.partner(current_atom)
                        else:
                            glyc_indices.append(bond.partner(current_atom).index) # add O bound to C2
                            CandO = getEsterCarbonyl(bond.partner(current_atom), current_atom)
                            carbonylC_indices.append(CandO[0])
                            glyc_indices.extend(CandO) # add C=O bound to C2O if present
                    glyc_indices.append(next_atom.index) # add third glycerol C
                    prev_atom = current_atom
                    current_atom = next_atom
                    if isGlycerolCarbon(2,current_atom):
                        glycerol_found = True
                        for bond in current_atom.bonds:
                            if bond.partner(current_atom).element == "O":
                                glyc_indices.append(bond.partner(current_atom).index) # add O bound to C3
                                CandO = getEsterCarbonyl(bond.partner(current_atom), current_atom)
                                carbonylC_indices.append(CandO[0])
                                glyc_indices.extend(CandO) # add C=O bound to C3O if present
        if not glycerol_found:
            glyc_indices = None
            carbonylC_indices = None
        return glyc_indices, carbonylC_indices


    phos_indices = hasPhosGroup(atomgroup)
    if phos_indices:
        also_has_glycerol_attached, ester_carbonyl_Cs = hasGlycerolAttached(atomgroup, phos_indices)
        if also_has_glycerol_attached and ester_carbonyl_Cs:
            return True
    return False


def main():
    """
    This script creates gmx workable groups and stores them in an ndx file.
    The groups are created based on whether or not they are phospholipids and reside above or below the center of geometry of all phospholipids.

    INPUT
        .tpr and .gro
    OUTPUT
        .ndx

    python3 memSplitNDX.py in.tpr in.gro -o index.ndx
    """
    args = parseArguments()
    running_folder = os.getcwd()

    if not os.path.exists(args.filename1):
        sys.exit("Error: " + "The file `" + running_folder +"/"+ args.filename1 + "` does not exist.")

    if not os.path.exists(args.filename2):
        sys.exit("Error: The file `" + running_folder + "/" + args.filename2 + "` does not exist.")

    # Create the Universe object
    u = mda.Universe(args.filename1, args.filename2)
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
    if len(phospholipids) == 0:
        sys.exit(f"No phospholipids in {args.filename2}.")

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
    if not args.force:
        for resname in phospholipids:
            assert above[resname] == below[resname]

    universe_groups = {} # {"<resname>": [<atomGroup>, <atomGroup>]}
    for resname in phospholipids:
        universe_groups[resname] = {}
        for z in z_plus_minus:
            str_resids = ' '.join([str(num) for num in phoslip_resids[resname][z]])
            universe_groups[resname][z] = u.select_atoms('resid ' + str_resids)

    # Output
    assert args.output.endswith(".ndx")
    if args.newNDX:
        writemode = 'w'
    elif os.path.exists(args.output):
        writemode = 'a'
    else:
        writemode = 'w'
    with mda.selections.gromacs.SelectionWriter(args.output, mode=writemode) as ndx:
        for resname in universe_groups:
            for z in z_plus_minus:
                ndx.write(universe_groups[resname][z], name=f"{resname}-" + z)
                ndx._outfile.write("\n")

    # Log
    print("\nCreated the following atom groups:")
    for resname in universe_groups:
        for groupname in universe_groups[resname]:
            print(f"    {resname}_{groupname}")
    print(f"\nThey were stored in:")
    print(f"    {args.output}")

main()
