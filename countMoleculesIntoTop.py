import warnings # MDAnalysis throws deprecationwarning on xdrlib which is slated for removal in Python 3.13
warnings.filterwarnings("ignore", message=".*'xdrlib' is deprecated.*", category=DeprecationWarning)

import xdrlib
import os
import sys
import argparse
import numpy
import MDAnalysis as mda
from collections import Counter
import re
from math import ceil

def parseArguments():
    parser = argparse.ArgumentParser(prog='countMoleculesIntoTop.py',
                                     description='Make the `[ molecules ]` block moleculenames in the second file (.top) match the residuenames present in the first file (.gro). Any moleculename not present in the second file is left unaltered.',
                                     epilog='Written by Lucas Roeleveld')

    parser.add_argument('coordfile',
                        help='the coordinate file')
    parser.add_argument('topfile',
                        help='the .top file')
    parser.add_argument('-o',
                        '--output',
                        default='topol.top',
                        help='specify a different output file')
    parser.add_argument('-sl',
                        '--skiplines',
                        default=0,
                        type=int,
                        help='set this to explicitly approve a number of lines inside the molecule block')
    parser.add_argument('-dontAdd',
                        action='store_true',
                        help='set this to refraim from counting molecules into the topology that arent there yet')
    parser.add_argument('-am',
                        '--addMols',
                        type=str,
                        default=None,
                        help='add lines to the molecules block with the given names')
    args = parser.parse_args()
    if args.addMols is not None:
        split_input = args.addMols.split()
        assert len(split_input) % 2 == 0 or len(split_input) == 1, "misformatted args.addMols argument, use either a single name to add a single molecule or space separated name-number pairs"
        if len(split_input) > 1:
            for i in range(0,len(split_input),2):
                try:
                    int(split_input[i+1])
                except:
                    sys.exit(f"`{args.addMols}` should be name-number pairs separated by spaces")
    return args


def main():
    """
    !!! args.addMols and addMols are missbehaving

    INPUT
        .gro and .top
    OUTPUT
        .top

    python3 countMoleculesIntoTop.py input.gro -p topol.top
    """
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASX", "GLX",
        "XAA", "UNK"
    ]
    spacer = 3

    args = parseArguments()
    if args.output is not None:
        outfile = args.output
    else:
        outfile = args.topfile

    running_folder = os.getcwd()

    if not os.path.exists(args.coordfile):
        sys.exit("Error: The file `" + running_folder + "/" + args.coordfile + "` does not exist.")
    if not os.path.exists(args.topfile):
        sys.exit("Error: The file `" + running_folder + "/" + args.topfile + "` does not exist.")

    if args.addMols is not None:
        split_input = args.addMols.split()
        addMols = {}
        if len(split_input) == 1:
            addMols[split_input[0]] = 1
        elif len(split_input) > 1:
            for i in range(0,len(split_input),2):
                addMols[split_input[i]] = int(split_input[i+1])

    # Create the Universe object
    u = mda.Universe(args.coordfile)
    mol_counter = Counter()
    mol_counter.update([residue.resname for residue in u.residues])
    for aa in amino_acids:
        if aa in mol_counter.keys():
            del mol_counter[aa]

    start_marker = r"\s*\[ molecules \]\s*"
    end_marker = r"\s*\[ (.*?) \]\s*"
    molecules_block = []  # Filled with the molecules block (also empty lines)
    mols_n_counts = [[], []]  # What is inside the .top file # one entry in both lists for each molecule line

    with open(args.topfile, "r") as file:
        lines = file.readlines()

    subset_indices = []  # indices to molecules_block
    # fill mols_n_counts
    inside_marker = 0
    j = 0
    # This fails when the `molecules` block is followed by a # instead of a [
    for i, line in enumerate(lines):
        if not inside_marker:
            if re.match(start_marker, line):
                molecules_block.append(line)
                inside_marker = i
        elif re.match(";", line) or re.match(r"\s*\n", line):
            j += 1
            molecules_block.append(line)
            continue
        elif re.match(end_marker, line):
            break
        else:
            j += 1
            molecules_block.append(line)
            subset_indices.append(j)
            mols_n_counts[0].append(line.split()[0])
            mols_n_counts[1].append(line.split()[1])


    # Check
    if len(mols_n_counts[0]) == 0 and len(mol_counter) == 0:
        sys.exit("There is nothing to count or change.")

    new_lines = []
    if j != 0:  # `[ molecules ]` block found
        # counting update
        assert args.skiplines in range(len(subset_indices) + 1), 'trying to skip an impossible number of lines in the molecules block.'
        remove_these_lines = []  # Will hold the line indices of the lines in the molcules block that should be discarded
        #j = 0
        if len(subset_indices) > 1:
            for i, index in enumerate(subset_indices):
                if i < args.skiplines: # Still need to empty mol_counter
                    if args.addMols is not None:
                        assert molecules_block[index].split()[0] not in addMols.keys(), f"You can not add `{molecules_block[index].split()[0]}` and also skip the line it is in"  # Is this line actually needed?
                    if molecules_block[index].split()[0] in mol_counter:
                        del mol_counter[molecules_block[index].split()[0]]
                    continue
                if mols_n_counts[0][i] in mol_counter:
                    if args.addMols and mols_n_counts[0][i] in addMols.keys():
                        mol_counter[mols_n_counts[0][i]] += addMols[mols_n_counts[0][i]]
                        del addMols[mols_n_counts[0][i]]
                    mols_n_counts[1][i] = mol_counter[mols_n_counts[0][i]]
                    del mol_counter[mols_n_counts[0][i]]
                else:
                    remove_these_lines.append(index)
            # remove the lines from the list
            for i in range(j + 1):
                del lines[inside_marker]


        # Strip empty lines from end of `molecules_block`
        for i in range(len(molecules_block)):
            if molecules_block[-1].lstrip(' ') == '\n':
                del molecules_block[-1]
            else:
                break

        ### Output formatting
        maxlen = [max([len(str_) for str_ in mols_n_counts[0]]) + spacer,
                  max([len(str(str_)) for str_ in mols_n_counts[1]])]
        for i, val in enumerate([14,9]):
            maxlen[i] = max(val, maxlen[i])
        for i, line in enumerate(molecules_block):
            if i in subset_indices:
                if maxlen[0] + maxlen[1] + 1 < len(line):
                    maxlen[0] = len(line) - maxlen[1] - 1
        if args.addMols:
            for key, val in addMols.items():
                if len(key) + spacer > maxlen[0]:
                    maxlen[0] = len(key) + spacer
                if len(str(val)) > maxlen[1]:
                    maxlen[1] = len(str(val))

        ###

        # make the new lines
        j = 0
        for i, line in enumerate(molecules_block):
            if i in remove_these_lines:  # Do nothing or comment the line
                if i in subset_indices:  # Need to increase j to access correct mols_n_counts
                    j += 1
                if args.topfile == outfile:
                    new_lines.append('; ' + line)
                    continue
            elif i in subset_indices:
                new_lines.append(mols_n_counts[0][j].ljust(maxlen[0], ' ') + str(mols_n_counts[1][j]).rjust(maxlen[1], ' ') + '\n')
                j += 1
            else:
                new_lines.append(line)
        if args.addMols:
            for key, val in addMols.items():
                new_lines.append(key.ljust(maxlen[0], ' ') + str(val).rjust(maxlen[1], ' ') + '\n')
            # Append other molecules counted in the file
            if not args.dontAdd and len(mol_counter) > 0:
                for key in mol_counter:
                    new_lines.append(key.ljust(maxlen[0], ' ') + str(mol_counter[key]).rjust(maxlen[1], ' ') + '\n')
        new_lines.append('\n')
    elif j == 0: # Creating a molecules block
        if args.skiplines:
            print(f"Warning: you tried to skip a line using the `-skiplines`/`-sl` flag, but there is no molecules block in {args.coordfile}")

        if args.addMols is not None:
            for key, val in addMols.items():
                if key in mol_counter.keys():
                    mol_counter[key] += val
                else:
                    mol_counter[key] = val
        maxlen = [max([len(str_) for str_ in mol_counter.keys()]) + spacer,
                  max([len(str(str_)) for str_ in mol_counter.values()])]
        for i, val in enumerate([14,9]):
            maxlen[i] = max(val, maxlen[i])
        if args.addMols is not None:
            for key, val in addMols.items():
                if len(key) + spacer > maxlen[0]:
                    maxlen[0] = len(key) + spacer
                if len(str(val)) > maxlen[1]:
                    maxlen[1] = len(str(val))
        new_lines.append('\n[ molecules ]\n')
        new_lines.append('; Compound'.ljust(maxlen[0], ' ') + ' #mols'.rjust(maxlen[1], ' ') + '\n')
            #new_lines.append(args.addMols.ljust(maxlen[0], ' ') + '1'.rjust(maxlen[1], ' ') + '\n')
        for key, count in mol_counter.items():
            new_lines.append(key.ljust(maxlen[0], ' ') + str(count).rjust(maxlen[1], ' ') + '\n')
    # Remove empty lines at end of file
        for i in reversed(range(len(lines))):
            if lines[i].strip(' ') == '\n':
                del lines[i]
            else:
                break
        inside_marker = len(lines)


    # Add new_lines to actual lines
    for i, line in enumerate(new_lines):
        lines.insert(inside_marker + i, line)

    # Remove empty lines at end of file
    for i in reversed(range(len(lines))):
        if lines[i].strip(' ') == '\n':
            del lines[i]
        else:
            break

    with open(outfile, 'w') as file:
        file.writelines(lines)


main()
