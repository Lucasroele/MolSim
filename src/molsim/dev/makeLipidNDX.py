###################
# makeLipidNDX.py #
################### 
#
# This program finds lipid residues and groups the non-hydrogen atoms in an mdanalysis universe which it reads from a .gro and a
# .tpr file into 5 groups per distinct lipid residue:
#     groupnames = ["phosphate", "glycerolester", "middle_chain", "outer_chain", "headgroup"]
# It then creates a .ndx file with all those groups.
#
# Slicing would look like this:
#
#     |   O   |     O  |
# HG -|- OPO -|- CCOC -|- middle_chain
#     |   O   |   |    |
#     |       |   COC -|- outer_chain
#     |       |     O  |
#
# This script should be ran using:
#     python3 makeLipidNDX.py in.gro in.tpr
#
# in which case it will generate
#     densitygroups_[<LIPID1>_<LIPID2>_...].ndx
#
# This program was written by Lucas Roeleveld
#


import os
import sys
import MDAnalysis as mda
import numpy
from collections import Counter
import argparse

# FUNCTIONS
# used to return the atomindex an atomindex is bound to given a bond between the two


def parseArguments():
    parser = argparse.ArgumentParser(prog='makeLipidNDX.py',
                                     description='Add lipid groups to existing .ndx or create new .ndx containing only lipid groups',
                                     epilog='Written by Lucas Roeleveld')

    parser.add_argument('filename1',
                        help='the coordinate file that contains lipids')
    parser.add_argument('filename2',
                        help='the .tpr file')
    parser.add_argument('-o',
                        '--output',
                        default='index.ndx',
                        help='the ndx file that will see the lipid groups added.')
    parser.add_argument('-uneven',
                        '--uneven',
                        action='store_true',
                        help='enable this flag when leaflet phospholipid distribution is not identical.')
    parser.add_argument('-wm',
                        '--writeMode',
                        type=str,
                        nargs='?',
                        help='can be set to write or append `w` or `a`, ovewrites any file already present when set to write.')
    parser.add_argument('-ll',
                        '--leaflet',
                        action='store_true',
                        help='Create all groups also for each leaflet')
    args = parser.parse_args()
    assert os.path.exists(args.filename1), f"Error: The file `{args.filename1}` does not exist."
    assert os.path.exists(args.filename2), f"Error: The file `{args.filename2}` does not exist."
    assert args.writeMode == 'a' and os.path.exists(args.output) or args.writeMode == 'w' or args.writeMode is None
    if not args.output.endswith(".ndx"):
        args.output += ".ndx"
    if args.writeMode is None:
        if os.path.exists(args.output):
            sys.exit(f"`{args.output}` exists, but you didn't specify the `-wm`.")
        else:
            args.writeMode = 'w'
    return args


def getOtherElement(bond, num): # Does the same as bond.partner(Atom) but on an index level
    if num == bond[0]:
        return bond[1]
    elif num == bond[1]:
        return bond[0]
    else:
        # Handle the case where the number is not found in the list
        raise ValueError("The provided number is not an element of the list.")

def makeGroupFromIndices(molmd, indices):
    str_indices = ' '.join([str(num) for num in indices])
    group = molmd.select_atoms('index ' + str_indices)
    return group

def findBondsToHydrogensInRes(u, resid):
    bonds_to_be_deleted = []
    molmd = u.select_atoms('resid ' + str(resid))
    for bond in molmd.bonds.indices:
        if molmd.select_atoms("index " + str(bond[0])).elements[0] == "H" or molmd.select_atoms("index " + str(bond[1])).elements[0] == "H":
            bonds_to_be_deleted.append(bond)
    return bonds_to_be_deleted

def hasPhosGroup(u, resid):
    has_phos = []
    molmd = u.select_atoms(f'resid {resid}')
    Ps = molmd.select_atoms('element P')
    if not Ps:
        return []
    else:
        P = Ps[0]
    Os = []
    for bond in P.bonds:
        if bond.partner(P).element == "O":
            Os.append(bond.partner(P).index)
    if len(Os) == 4:
        has_phos.append(P.index)
        has_phos.extend(Os)
    return has_phos

# Check if atomindex in molmd is a O-C-C or O-C-2C
def isGlycerolCarbon(carbontype, atom):
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

def getEsterCarbonyl(bridgeO, alcoholC):
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
                    return []
                else:
                    other_atoms += 1
        else:
            return []
    return carbonyl_indices

def hasGlycerolAttached(u, phos_indices):
    glycerol_found = False
    carbonylC_indices = []
    current_atom_index = 0 # the index of the atom the decider is at
    phosmol = makeGroupFromIndices(u, phos_indices)
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
        glyc_indices = []
    return glyc_indices, carbonylC_indices

# The parseSegment fuction is not particularyly good code and should probably be rewrtitten
# to use bond.partner(Atom) and not user_defined function getOtherElement

# traverse through molecule and assign atoms to groups
# Fetches all atoms adjacent and connected to current_index not in grouped_indices
def parseSegment(current_index, molmd, grouped_indices):
    more_molecules = False
    grouped_indices_array = []
    for key, value in grouped_indices.items():
        grouped_indices_array.extend(value)
    visited = []
    forks = {} # {'<index_that_forks>': [<indices_it_forks_to>]}
    isfork = [] # when len(isfork) >= 2 -> current_index is fork
    for bond in molmd.select_atoms("index " + str(current_index)).bonds.indices:
        index_under_cons = getOtherElement(bond, current_index)
        if index_under_cons in grouped_indices_array:
            continue
        else:
            isfork.append(index_under_cons)
            more_molecules = True
    if len(isfork) > 1:
        forks[str(current_index)] = isfork[0:-1]
        next_index = isfork[-1]
    elif len(isfork) == 1:
        next_index = isfork[0]

    while more_molecules:
        isfork = []
        current_index = next_index
        visited.append(current_index)

        # find indices connected to current_index
        for bond in molmd.select_atoms("index " + str(current_index)).bonds.indices:
            index_under_cons = getOtherElement(bond, current_index)
            if index_under_cons in grouped_indices_array or index_under_cons in visited:
                continue
            else:
                isfork.append(index_under_cons)

        # Decide where to go
        more_molecules = True
        if len(isfork) >= 2: # add fork to forkdict
            forks[str(current_index)] = isfork[0:-1]
            next_index = isfork[-1]
        elif len(isfork) == 1: # just continue along chain
            next_index = isfork[0]
        elif len(forks) >= 1: # go back to last fork
            popped_item = forks.popitem()
            if len(popped_item[1]) > 1:
                forks[popped_item[0]] = popped_item[1][0:-1]
            next_index = popped_item[1][-1]
        else:
            more_molecules = False
    return visited

args = parseArguments()

# Create the Universe object
u = mda.Universe(args.filename2, args.filename1)

# Count number of each residue present
residue_counter = Counter()
residue_counter.update([residue.resname for residue in u.residues])
resnames = [key for key, value in residue_counter.items()]

# Log
print(f"Found {len(residue_counter)} residues:")
for key in residue_counter:
    print(f"    {key}")

# Make `phospholipids` hold the resnames to study
resnames_with_phosphate = [] # [""]
has_phos = {} # {"<resname>": [#]}
also_has_glycerol_attached = {} # {"<resname>": [#]}
ester_carbonyl_Cs = {} # {"<resname>": [#]}
phospholipids = [] # ["<resname>"]
resids = {} # {"<resname>": "<resid>"}

# Figure out which resnames belong to phospholipids
for resname in resnames:
    for residue in u.residues:
        if residue.resname == resname:
            resids[resname] = residue.resid
            break
    has_phos[resname] = hasPhosGroup(u, resids[resname])
    if has_phos[resname]:
        resnames_with_phosphate.append(resname)
for resname in resnames_with_phosphate:
    also_has_glycerol_attached[resname], ester_carbonyl_Cs[resname] = hasGlycerolAttached(u, has_phos[resname])
    if also_has_glycerol_attached[resname]:
        phospholipids.append(resname)
if len(phospholipids) == 0:
    sys.exit("\nCould not find any phospholipids. This can happen when the phospholipid contains multiple P atoms.")
else:
    print(f"\n{len(phospholipids)} of these were identified as phospholipids:")
    for phospholipid in phospholipids:
        print(f"    {phospholipid}")

groupnames = ["phosphate", "glycerolester", "middle_chain", "outer_chain", "headgroup"]
molmd = {} # {"<resname>": <atomGroup>}
resname_groups = {} # {"resname": {"groupname": <atomGroup>}}
grouped_indices = {} # {"resname": [#]}
all_grouped_indices = {} # {"resname": [#]}
there_was_a_problem = {} # {"<resname>": <Bool>}
universe_group_indices = {} # {"<resname>": {"groupname": [#]}}
universe_groups = {} # {"<resname>": {"groupname": <atomGroup>}}
for resname in phospholipids:
    resname_groups[resname] = {}
    grouped_indices[resname] = {}
    bonds_to_be_deleted = findBondsToHydrogensInRes(u, resids[resname])
    u.delete_bonds(bonds_to_be_deleted)
    molmd[resname] = u.select_atoms('resid ' + str(resids[resname]) + " and not element H")
    resname_groups[resname][groupnames[0]] = makeGroupFromIndices(molmd[resname], has_phos[resname])
    grouped_indices[resname][groupnames[0]] = has_phos[resname]
    resname_groups[resname][groupnames[1]] = makeGroupFromIndices(molmd[resname], also_has_glycerol_attached[resname])
    grouped_indices[resname][groupnames[1]] = also_has_glycerol_attached[resname]
    grouped_indices[resname][groupnames[2]] = parseSegment(ester_carbonyl_Cs[resname][0], molmd[resname], grouped_indices[resname])
    resname_groups[resname][groupnames[2]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[2]])
    grouped_indices[resname][groupnames[3]] = parseSegment(ester_carbonyl_Cs[resname][1], molmd[resname], grouped_indices[resname])
    resname_groups[resname][groupnames[3]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[3]])
    for atom in resname_groups[resname][groupnames[0]]:
        if atom.element == "P":
            continue
        elif len(atom.bonds) == 1:
            continue
        else:
            for bond in atom.bonds:
                if not bond.partner(atom) in resname_groups[resname][groupnames[0]] and not bond.partner(atom) in resname_groups[resname][groupnames[1]]:
                    grouped_indices[resname][groupnames[4]] = parseSegment(atom.index, molmd[resname], grouped_indices[resname])
                    resname_groups[resname][groupnames[4]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[4]])
                    break

    all_grouped_indices[resname] = []
    for key in grouped_indices[resname]:
        all_grouped_indices[resname] += grouped_indices[resname][key]
    errmes = f"Something went wrong in grouping the backbone atoms of {resname}."
    assert len(all_grouped_indices[resname]) == len(molmd[resname].indices), errmes
    assert len(set(molmd[resname].indices) - set(all_grouped_indices[resname])) == 0, errmes
    assert len(all_grouped_indices[resname]) == len(set(all_grouped_indices[resname])), errmes

#    if len(all_grouped_indices[resname]) != len(molmd[resname].indices):
#        there_was_a_problem[resname] = True
#    elif len(set(molmd[resname].indices) - set(all_grouped_indices[resname])) > 0:
#        there_was_a_problem[resname] = True
#    elif len(all_grouped_indices[resname]) - len(set(all_grouped_indices[resname])) > 0:
#        there_was_a_problem[resname] = True
#    else:
#        there_was_a_problem[resname] = False
#
#    if there_was_a_problem[resname]:
#        sys.exit(errmes)


    # Group all resnames in the universe
    universe_group_indices[resname] = {}
    for groupname in resname_groups[resname].keys():
        universe_group_indices[resname][groupname] = []
    for atom in u.atoms:
        if atom.element == "H":
            continue
        for groupname, group in resname_groups[resname].items():
            if atom.resname == resname and atom.name in group.names:
                universe_group_indices[resname][groupname].append(atom.index)
                break
    universe_groups[resname] = {}
    for groupname in resname_groups[resname].keys():
        universe_groups[resname][groupname] = makeGroupFromIndices(u, universe_group_indices[resname][groupname])
    for groupname in universe_group_indices[resname]:
        assert len(universe_group_indices[resname][groupname])/residue_counter[resname] == len(grouped_indices[resname][groupname]), errmes
#        if len(universe_group_indices[resname][groupname])/residue_counter[resname] != len(grouped_indices[resname][groupname]):
#            there_was_a_problem[resname] = True
#
#    if there_was_a_problem[resname]:
#        sys.exit(f"Something went wrong in grouping the backbone atoms of {resname}.")

universe_groups['lipid'] = makeGroupFromIndices(u, [c for a in universe_group_indices.values() for b in a.values() for c in b])

# Output
with mda.selections.gromacs.SelectionWriter(args.output, mode=args.writeMode) as ndx:
    for resname in universe_groups:
        if resname == 'lipid':
            ndx.write(universe_groups[resname], name=f"{resname}")
            ndx._outfile.write("\n")
        else:
            for groupname, group in universe_groups[resname].items():
                ndx.write(group, name=f"{resname}_{groupname}")
                ndx._outfile.write("\n")

# Log
print("\nCreated the following atom groups:")
for resname in universe_groups:
    if resname == 'lipid':
        print(f"    {resname}")
    else:
        for groupname in universe_groups[resname]:
            print(f"    {resname}_{groupname}")
print(f"\nThey were stored in:")
print(f"    {args.output}")


