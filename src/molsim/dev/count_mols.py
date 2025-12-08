import sys
import argparse
import MDAnalysis as mda
from collections import Counter
import re
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('count_mols',
                                   help='Make the `[ molecules ]` block moleculenames in the second file (.top) match the residuenames present in the first file (.gro). Any moleculename not present in the second file is left unaltered')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='count_mols.py',
                                     description='Make the `[ molecules ]` block moleculenames in the second file (.top) match the residuenames present in the first file (.gro). Any moleculename not present in the second file is left unaltered.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('coordfile',
                        help='the coordinate file')
    parser.add_argument('topfile',
                        help='the .top file')
    parser.add_argument('-o',
                        '--output',
                        help='specify a different output file')
    parser.add_argument('-sl',
                        '--skiplines',
                        default=None,
                        type=int,
                        help='Set this to explicitly approve a number of lines inside the molecule block. Useful for macromolecules.')
    parser.add_argument('-dontAdd',
                        action='store_true',
                        help='Set this to refraim from counting molecules into the topology that arent there yet.')
    parser.add_argument('-am',
                        '--addMols',
                        type=str,
                        default=None,
                        help='Add extra molecules to the molecules block with the given names.')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='Overwrite the output file if it already exists.')
    return parser


def validateArguments(args):
    if args.addMols is not None:
        split_input = args.addMols.split()
        len_split = len(split_input)
        assert len_split % 2 == 0 or len_split == 1, "misformatted args.addMols argument, use either a single name to add a single molecule or space separated name-number pairs"
        if len_split > 1:
            for i in range(0,len_split,2):
                assert split_input[i+1].isdigit(), f"`{args.addMols}` should be a string of name-number pairs separated by spaces"
    assert Path(args.coordfile).is_file(), f"Coordinate file `{args.coordfile}` does not exist"
    assert Path(args.topfile).is_file(), f"Topology file `{args.topfile}` does not exist"
    if args.output is None:
        args.output = args.topfile
    else:
        assert not Path(args.output).is_file() or args.force, f"Output file `{args.output}` already exists and `-f` was not used"
    return args


def main(args):
    """
    !!! args.addMols and addMols are missbehaving

    INPUT
        .gro and .top
    OUTPUT
        .top

    python3 count_mols.py input.gro topol.top
    """
    args = validateArguments(args)
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASX", "GLX",
        "XAA", "UNK"
    ]
    spacer = 3


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
    # old end_marker = r"\s*\[ (.*?) \]\s*"
    end_marker = r"\s*\[ (.*?) \]\s*|\s*#.*"
    molecules_block = []  # Filled with the molecules block (also empty lines)
    mols_n_counts = [[], []]  # What is inside the .top file, one entry in both lists for each molecule line

    with open(args.topfile, "r") as file:
        lines = file.readlines()

    subset_indices = []  # indices to lines inside molecules_block that have data
    # fill mols_n_counts
    inside_marker = 0
    comment_j = []
    j = 0
    for i, line in enumerate(lines):
        if not inside_marker:
            if re.match(start_marker, line):
                molecules_block.append(line)
                inside_marker = i
        elif re.match(";", line) or re.match(r"\s*\n", line):
            j += 1
            molecules_block.append(line)
            comment_j.append(j)
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
    if inside_marker:  # `[ molecules ]` block found
        # counting update
        assert args.skiplines is None or args.skiplines in range(len(subset_indices) + 1), 'Trying to skip more lines than there are in the `[ molecules ]` block'
        remove_these_lines = []  # Will hold the line indices of the lines in the molcules block that should be discarded
        if len(subset_indices) > 1:
            for i, index in enumerate(subset_indices):
                if args.skiplines is not None and i < args.skiplines: # Still need to empty mol_counter
                    if args.addMols is not None:
                        assert molecules_block[index].split()[0] not in addMols.keys(), f"You can not add `{molecules_block[index].split()[0]}` and also skip the line it is in"  # Is this line actually needed?
                    if molecules_block[index].split()[0] in mol_counter:
                        del mol_counter[molecules_block[index].split()[0]]
                    continue
                if mols_n_counts[0][i] in mol_counter:
                    # adding to existing lines
                    if args.addMols is not None and mols_n_counts[0][i] in addMols.keys():
                        mol_counter[mols_n_counts[0][i]] += addMols[mols_n_counts[0][i]]
                        del addMols[mols_n_counts[0][i]]
                    mols_n_counts[1][i] = mol_counter[mols_n_counts[0][i]]
                    del mol_counter[mols_n_counts[0][i]]
                else: # the line present inside `*.top` is bs and should be removed
                    remove_these_lines.append(index)
            # remove the `[ molecules ]` block lines inside the .top file
            for i in range(j + 1):
                del lines[inside_marker]


        # Strip empty lines from end of `molecules_block`
        for i in range(len(molecules_block)):
            if molecules_block[-1].lstrip(' ') == '\n':
                del molecules_block[-1]
            else:
                break

        ### Output formatting
        ## entirely for reading header comment which default to '; Compound #mols'
        header_comment_j = None
        if comment_j is not None:
            for j in comment_j:
                if min(subset_indices) - j == 1:
                    if len(molecules_block[j].split()) == 3:
                        header_comment_j = j
                        header_comment = molecules_block[j].rstrip('\n')
                        header_comment_len = [len(' '.join(header_comment.split()[0:1])) + spacer, len(header_comment.split()[-1])]
        
        ## assessing justification
        maxlen = [max([len(str_) for str_ in mols_n_counts[0]]) + spacer,
                  max([len(str(str_)) for str_ in mols_n_counts[1]])]
        for i, val in enumerate([10,5]): # 10 for molecule name, 5 for number and '#mols'
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
        if header_comment_j is not None:
            for i, val in enumerate(header_comment_len):
                if val > maxlen[i]:
                    maxlen[i] = val
        ##
        ## Creating the header comment
        if header_comment_j is not None:
            header_comment = (header_comment.split()[0] + ' ' + \
                              header_comment.split()[1]).ljust(maxlen[0], ' ') + \
                             header_comment.split()[2].rjust(maxlen[1], ' ') + '\n'
        else:
            header_comment = '; Compound'.ljust(maxlen[0], ' ') + '#mols'.rjust(maxlen[1], ' ') + '\n'
        ##
        ###

        # make the new lines
        j = 0
        for i, line in enumerate(molecules_block):
            if i in remove_these_lines:  # Do nothing or comment the line
                if i in subset_indices:  # Need to increase j to access correct mols_n_counts
                    j += 1
                if args.topfile == args.output: # might actually need the line
                    new_lines.append('; ' + line)
            elif i in subset_indices:
                if i == min(subset_indices) and header_comment_j is None:
                    new_lines.append(header_comment)
                new_lines.append(mols_n_counts[0][j].ljust(maxlen[0], ' ') + str(mols_n_counts[1][j]).rjust(maxlen[1], ' ') + '\n')
                j += 1
            elif header_comment_j is not None and i == header_comment_j:
                new_lines.append(header_comment)
            else:
                new_lines.append(line)
        if args.addMols:
            for key, val in addMols.items():
                if key in mol_counter.keys():
                    mol_counter[key] += val
                    continue
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

    with open(args.output, 'w') as file:
        file.writelines(lines)

if __name__ == '__main__':
    args = parseArguments()
    main(args)
