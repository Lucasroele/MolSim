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


from pathlib import Path
import MDAnalysis as mda
from collections import Counter
import argparse
import itertools as it
#from molsim.utils.mda import isPhospholipid
#from molsim.parsers import NdxParser
from molsim.utils.mda import hasPhosGroup
#from molsim.utils.mda import getEsterCarbonyl
#from molsim.utils.mda import isGlycerolCarbon
from molsim.utils.mda import hasGlycerolAttached
from molsim.utils.mda import leafletResolution

# FUNCTIONS
def register(subparsers):
    parser = subparsers.add_parser('lipid_ndx',
                                   help='Add lipid groups to existing .ndx or create new .ndx containing only lipid groups')
    parser = addArguments(parser)
    parser.set_defaults(func=main)

        
def parseArguments():
    parser = argparse.ArgumentParser(prog='lipid_ndx.py',
                                     description='Add lipid groups to existing .ndx or create new .ndx containing only lipid groups',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('topfile',
                        help='the .tpr file')
    parser.add_argument('-o',
                        '--output',
                        default='index.ndx',
                        help='the ndx file that will see the lipid groups added.')
    parser.add_argument('-u',
                        '--uneven',
                        action='store_true',
                        help='enable this flag when leaflet phospholipid distribution is not identical.')
    parser.add_argument('-a',
                        '--append',
                        action="store_true",
                        help='set to append the groups to the ndx file.')
    parser.add_argument('-ll',
                        '--leaflet',
                        action='store_true',
                        help='Create all groups also for each leaflet')
    parser.add_argument('-nos',
                        '--no_subgroups',
                        action='store_true',
                        help='disable subgroups, i.e. phosphate, glycerol and glycerolester')
    parser.add_argument('-q',
                        '--quiet',
                        action='store_true',
                        help='disable logging')
    return parser


def validateArguments(args):
    if not args.output.endswith(".ndx"):
        args.output += ".ndx"
    outpath = Path(args.output)
    toppath = Path(args.topfile)
    assert toppath.is_file(), f"Error: The file `{args.topfile}` does not exist."
    if args.append:
        assert outpath.is_file(), f"Error: The file `{args.output}` does not exist."
    return args

def makeGroupFromIndices(molmd, indices):
    str_indices = ' '.join([str(num) for num in indices])
    group = molmd.select_atoms('index ' + str_indices)
    return group


def findBondsToHydrogensInRes(atomgroup):
    bonds_to_be_deleted = []
    for bond in atomgroup.bonds.indices:
        if atomgroup.select_atoms("index " + str(bond[0])).elements[0] == "H" or atomgroup.select_atoms("index " + str(bond[1])).elements[0] == "H":
            bonds_to_be_deleted.append(bond)
    return bonds_to_be_deleted

# The parseSegment fuction is not particularyly good code and should probably be rewrtitten

# traverse through molecule and assign atoms to groups
# Fetches all atoms adjacent and connected to current_index not in grouped_indices

def parseSegment(current_index, atomgroup, grouped_indices: list):
    """
    traverse through molecule and assign atoms to groups
    

    Parameters
    ----------
    current_index :
        should be not inside grouped_indices and will be part of the returned group
    atomgroup :
        atomgroup
    grouped_indices :
        grouped_indices

    Returns
    -------
    None.
    """
    more_atoms = False
    visited = []
    forks = {} # {'<index_that_forks>': [<indices_it_forks_to>]}
    isfork = [] # when len(isfork) >= 2 -> current_index is fork
    for bond in atomgroup.select_atoms("index " + str(current_index)).bonds:
        index_under_cons = bond.partner(atomgroup.universe.atoms[current_index]).index
        if index_under_cons in grouped_indices:
            continue
        else:
            isfork.append(index_under_cons)
            more_atoms = True
    if len(isfork) > 1:
        forks[str(current_index)] = isfork[0:-1]
        next_index = isfork[-1]
    elif len(isfork) == 1:
        next_index = isfork[0]

    while more_atoms:
        isfork = []
        current_index = next_index
        visited.append(current_index)

        # find indices connected to current_index
        for bond in atomgroup.select_atoms("index " + str(current_index)).bonds:
            index_under_cons = bond.partner(atomgroup.universe.atoms[current_index]).index
            if index_under_cons in grouped_indices or index_under_cons in visited:
                continue
            else:
                isfork.append(index_under_cons)

        # Decide where to go
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
            more_atoms = False
    return visited


def main(args):
    args = validateArguments(args)
    # Create the Universe object
    u = mda.Universe(args.topfile)

    # Count number of each residue present
    residue_counter = Counter()
    residue_counter.update([residue.resname for residue in u.residues])
    resnames = [key for key, value in residue_counter.items()]

    # Log
    if not args.quiet:
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
        # select a template phospholipid to work with
        for residue in u.residues:
            if residue.resname == resname:
                resids[resname] = residue.resid
                break
        molecule = u.select_atoms('resid ' + str(resids[resname]))
        has_phos[resname] = hasPhosGroup(molecule)
        if has_phos[resname] is not None:
            resnames_with_phosphate.append(resname)
    for resname in resnames_with_phosphate:
        also_has_glycerol_attached[resname], ester_carbonyl_Cs[resname] = hasGlycerolAttached(u, has_phos[resname])
        if also_has_glycerol_attached[resname]:
            phospholipids.append(resname)
    if len(phospholipids) == 0:
        if not args.quiet:
            print("\nCould not find any phospholipids. This can happen for example when the phospholipid contains multiple P atoms.")
        return
    else:
        if not args.quiet:
            print(f"\n{len(phospholipids)} of these were identified as phospholipids:")
            for phospholipid in phospholipids:
                print(f"    {phospholipid}")

    # Parse each one instance of each phospholipid
    groupnames = ["", "phosphate", "glycerolester", "middle_chain", "outer_chain", "headgroup"]
    molmd = {} # {"<resname>": <atomGroup>}
    resname_groups = {} # {"resname": {"groupname": <atomGroup>}}
    grouped_indices = {} # {"resname": {"groupname": [#]}}
    all_grouped_indices = {} # {"resname": [#]}
    universe_group_indices = {} # {"<resname>": {"groupname": [#]}}
    universe_groups = {} # {"<resname>": {"groupname": <atomGroup>}}
    for resname in phospholipids:
        resname_groups[resname] = {}
        grouped_indices[resname] = {}
        bonds_to_be_deleted = findBondsToHydrogensInRes(u.select_atoms('resid ' + str(resids[resname])))
        u.delete_bonds(bonds_to_be_deleted)
        molmd[resname] = u.select_atoms('resid ' + str(resids[resname]) + " and not element H")
        # phosphate
        resname_groups[resname][groupnames[1]] = makeGroupFromIndices(molmd[resname], has_phos[resname])
        grouped_indices[resname][groupnames[1]] = has_phos[resname]
        # glycerolester
        resname_groups[resname][groupnames[2]] = makeGroupFromIndices(molmd[resname], also_has_glycerol_attached[resname])
        grouped_indices[resname][groupnames[2]] = also_has_glycerol_attached[resname]
        # middle_chain
        grouped_indices[resname][groupnames[3]] = parseSegment(ester_carbonyl_Cs[resname][0], molmd[resname], list(it.chain.from_iterable(grouped_indices[resname].values())))
        resname_groups[resname][groupnames[3]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[3]])
        # outer_chain
        grouped_indices[resname][groupnames[4]] = parseSegment(ester_carbonyl_Cs[resname][1], molmd[resname], list(it.chain.from_iterable(grouped_indices[resname].values())))
        resname_groups[resname][groupnames[4]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[4]])
        # headgroup
        for atom in resname_groups[resname][groupnames[1]]:
            if atom.element == "P":
                continue
            elif len(atom.bonds) == 1:
                continue
            else:
                for bond in atom.bonds:
                    if bond.partner(atom) not in resname_groups[resname][groupnames[1]] and bond.partner(atom) not in resname_groups[resname][groupnames[2]]:
                        grouped_indices[resname][groupnames[5]] = parseSegment(atom.index, molmd[resname],list(it.chain.from_iterable(grouped_indices[resname].values())))
                        resname_groups[resname][groupnames[5]] = makeGroupFromIndices(molmd[resname], grouped_indices[resname][groupnames[5]])
                        break
        # phoslip
        resname_groups[resname][groupnames[0]] = molmd[resname]
        grouped_indices[resname][groupnames[0]] = list(molmd[resname].indices)

        # Sanity checks
        all_grouped_indices[resname] = []
        for key in grouped_indices[resname]:
            all_grouped_indices[resname].extend(grouped_indices[resname][key])
        errmes = f"Something went wrong in grouping the non-H atoms of {resname}."
        assert len(all_grouped_indices[resname]) == 2 * len(molmd[resname].indices), errmes + ' ' + str(len(all_grouped_indices[resname])) + ' vs ' + str(len(molmd[resname].indices))
        assert len(set(molmd[resname].indices) - set(all_grouped_indices[resname])) == 0, errmes
        assert len(all_grouped_indices[resname]) == 2 * len(set(all_grouped_indices[resname])), errmes

        # Group all atoms inside the universe not just the template residue atoms
        universe_group_indices[resname] = {}
        for groupname in resname_groups[resname].keys():
            universe_group_indices[resname][groupname] = []
        for atom in u.select_atoms('not element H').atoms:
            for groupname, group in resname_groups[resname].items():
                if atom.resname == resname and atom.name in group.names:
                    universe_group_indices[resname][groupname].append(atom.index)
        universe_groups[resname] = {}
        for groupname in resname_groups[resname].keys():
            universe_groups[resname][groupname] = makeGroupFromIndices(u, universe_group_indices[resname][groupname])
        for groupname in universe_group_indices[resname]:
            assert len(universe_group_indices[resname][groupname])/residue_counter[resname] == len(grouped_indices[resname][groupname]), errmes

    universe_groups['lipid'] = makeGroupFromIndices(u, [c for a in universe_group_indices.values() for b in a.values() for c in b])
    
    if args.no_subgroups:
        for resname in phospholipids:
            for groupname in groupnames[1:]:
                del universe_groups[resname][groupname]

    def nameMaker(nested_group_dict):
        """
        Recursive function that formats the nested groups

        Parameters
        ----------
        nested_group_dict :
            e.g. {"POPC": {"phophate": atomgroup1, "glycerol": atomgroup2, ...}, ...}

        Returns
        -------
        names :
            e.g. ["POPC_phosphate", "POPC_glycerol", ...]
        atoms :
            e.g. [atomgroup1, atomgroup2, ...]
        """
        names = []
        atoms = []
        for key1, val1 in nested_group_dict.items():
            if isinstance(val1, dict):
                temp_names, temp_atoms = nameMaker(val1)
                for key2, val2 in zip(temp_names, temp_atoms):
                    if key2 == "":
                        names.append(key1)
                    else:
                        names.append(f"{key1}_{key2}")
                    atoms.append(val2)
            else:
                names.append(key1)
                atoms.append(val1)
        return names, atoms

    names, atoms = nameMaker(universe_groups)
    if args.leaflet:
        new_names = []
        new_atoms = []
        for i in range(len(names)):
            leaflet_p, leaflet_n = leafletResolution(atoms[i])
            new_names.append(f"{names[i]}_Z+")
            new_names.append(f"{names[i]}_Z-")
            if len(leaflet_p) != len(leaflet_n) and not args.uneven:
                raise ValueError(f"Leaflet resolution failed for {names[i]}, set `--uneven` to ignore this error.")
            new_atoms.append(leaflet_p)
            new_atoms.append(leaflet_n)
        names = new_names
        atoms = new_atoms

    
    # Add the checking part for already present from split_mem

    # Output
    if args.append:
        write_mode = 'a'
        # Make sure the file ends with a empty line before appending
        with open(args.output, 'r+') as ndx:
            if ndx.readlines()[-1] != '\n':
                ndx.write("\n")
    else:
        write_mode = 'w'
    with mda.selections.gromacs.SelectionWriter(args.output, mode=write_mode) as ndx:
        for i, name in enumerate(names):
            ndx.write(atoms[i], name=name)
            ndx._outfile.write("\n")

    # Log
    if not args.quiet:
        print("\nCreated the following atom groups:")
        for name in names:
            print(f"    {name}")
        print("\nThey were stored in:")
        print(f"    {args.output}")

if __name__ == "__main__":
    args = parseArguments()
    main(args)
