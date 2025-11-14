import sys
import os
import argparse
from collections import OrderedDict

#
# This script only works when the topology file has:

# Requires terminal to be
#   N
#   H
#   C
#   ...
#   C
#   O

# Requires start to be
#   N
#   HN
#   CA
#   HCA
#   _CA
#   ...
#   C
#   O
def main():
    def parseArguments():
        default_out = 'topol.top'
        parser = argparse.ArgumentParser(prog='topologyCyclicizer.py',
                                         description='This script makes a peptide chain in a topology a cyclic peptide.',
                                         epilog='This script was written by Daria de Raffele and modified by Lucas Roeleveld.')

        parser.add_argument('inputFile',
                            help='The topology file of a peptide')           # positional argumen 1
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            default=default_out,
                            nargs='?',
                            help=f"name of the output which defaults to `{default_out}`")
        parser.add_argument("-f",
                            "--force",
                            action='store_true',
                            help='overrides some built-in safeguards')
        args = parser.parse_args()
        assert os.path.exists(args.inputFile), f"`{args.inputFile}` doesn't exist."
        if not args.force:
            assert args.inputFile.endswith(".top"), "`{args.inputFile}` does not have the .top extension; use `-f` to force"
        if os.path.exists(args.output) and not args.force:
            sys.exit(f"The outputfile `{args.output}` already exists, use `-f` to overwrite.")
        return args


    args = parseArguments()

    f = open(args.inputFile, 'r')

    letter = ['N', 'CA', 'C']

    l = f.readlines()   # lines
    b = []              # bond lines header line #
    w = []              # dihedral lines header line #
    q = []              # angle lines header line #
    z = []              # cmap lines header line #

    # row numbers where these strings are reported
    for dih, j in enumerate(l):
        if '[ bonds ]' in j:
            b.append(dih)
        if '[ dihedrals ]' in j:
            w.append(dih)
        if '[ angles ]' in j:
            q.append(dih)
        if '[ cmap ]' in j:
            z.append(dih)

    # generate the list of the peptide residues
    topo = []  # atomnum, resnum, atomtype
    backbone = []  # [[<atomnum>, <atomtype>]]
    N = []  # [<atomnum>]
    CA = []  # [<atomnum>]
    C = []  # [<atomnum>]
    i = 0
    while '[ bonds ]' not in l[i] and i < len(l):
        d = l[i].split()
        if not l[i].split():
            i += 1
            continue
        if d[0] == ';' or d[0] == '[':
            i += 1
            continue
        if d[0].isdigit():
            res = d[0], d[2], d[4]  # atomnum, resnum, atomtype
            topo.append(res)
        if d[0].isdigit() and d[4] in letter:
            bb = d[0], d[4]
            backbone.append(bb)
        i += 1
    for i in range(0, len(backbone)):  # CMAP!!!
        if backbone[i][1] == 'C':
            C.append(backbone[i][0])
        if backbone[i][1] == 'CA':
            CA.append(backbone[i][0])
        if backbone[i][1] == 'N':
            N.append(backbone[i][0])

    # This is added to move #include blocks and general information blocks at the bottom of the file to the top
    # It is assumed that all block of information are empty line separated
    # The headers specified in the next line remain at the same place in the file
    topolblocks = ['[ atoms ]', '[ bonds ]', '[ angles ]', '[ pairs ]', '[ dihedrals ]', '[ impropers ]', '[ cmap ]', '[ moleculetype ]']
    paragraphs = OrderedDict()
    par_order_in = []
    prologue_end = 0
    i = 0
    while i < len(l):
        prev_line_empty = False
        if 'prologue' not in par_order_in:
            if i < len(l) and not l[i].strip(' ').startswith('['):
                par_order_in.append('prologue')
                paragraphs['prologue'] = []
            while i < len(l) and not l[i].strip(' ').startswith('['):
                paragraphs['prologue'].append(l[i])
                i += 1
            prologue_end = i
        # Normal blocks
        if i < len(l) and l[i].strip(' ').startswith('['):
            block = l[i].rstrip('\n').strip(' ')
            if block in par_order_in:
                block += '_'  # Allows for blocks with the same name
            par_order_in.append(block)
            paragraphs[block] = []
            while i < len(l) and l[i].rstrip('\n').strip(' ') != '':
                paragraphs[block].append(l[i])
                i += 1
            i += 1
        # # and ; blocks
        elif i < len(l) and (l[i].strip(' ').startswith('#') or l[i].strip(' ').startswith(';')):
            i_start = i
            while i < len(l):
                if l[i].strip(' ').startswith('#ifdef') or l[i].strip(' ').startswith('#include'):
                    block = l[i].rstrip('\n').strip(' ')
                    i = i_start
                    break
                elif l[i].rstrip('\n').strip(' ') == '' or i + 1 == len(l):
                    block = l[i_start].rstrip('\n').strip(' ')
                    i = i_start
                    break
                i += 1
            par_order_in.append(block)
            paragraphs[block] = []
            while i < len(l) and l[i].rstrip('\n').strip(' ') != '':
                paragraphs[block].append(l[i])
                i += 1
            i += 1


    remove_these = []
    # remove topolblocks, they will just be read from the input file
    for i, parname in enumerate(par_order_in):
        for topolblock in topolblocks:
            if parname.startswith(topolblock):
                del paragraphs[parname]
                remove_these.append(i)
                break
    for index in reversed(remove_these):
        del par_order_in[index]

    # generate list of the first and last residues
    first = topo[0][1]  # num first resid
    last = topo[-1][1]  # num last resid
    fres = []  # atoms in fiirst resid
    lres = []  # atoms in last resid
    for i in range(len(topo)):
        if topo[i][1] == first:
            fres.append(topo[i])
        if topo[i][1] == last:
            lres.append(topo[i])

    assert fres[0][2][0] == 'N' and fres[1][2][0] == 'H' and fres[2][2][0] == 'C' and fres[-2][2][0] == 'C', "N-terminal should be N-H-C-...-C-O"
    assert lres[-1][2][0] == 'O' and lres[-2][2][0] == 'C' and lres[0][2][0] == 'N' and lres[1][2][0] == 'H' and lres[2][2][0] == 'C', "C-terminal should be N-H-C-...-C-O"

    # Lines that will be added
    # bond
    # Add bond between C being the penultimate atom in lres and N being the first in fres
    bond = '%5s %5s %12s' % (lres[-2][0], fres[0][0], '1 ;cycle')

    # pairs
    pair1 = '%5s %5s %12s' % (lres[0][0], fres[0][0], '1 ;cycle')  # N-N
    pair2 = '%5s %5s %12s' % (lres[2][0], fres[1][0], '1 ;cycle')  # CA-HN
    pair3 = '%5s %5s %12s' % (lres[2][0], fres[2][0], '1 ;cycle')  # CA-CA
    pair4 = '%5s %5s %12s' % (lres[3][0], fres[0][0], '1 ;cycle')  # HCA-N
    pair5 = '%5s %5s %12s' % (lres[4][0], fres[0][0], '1 ;cycle')  # _CA-N
    pair6 = '%5s %5s %12s' % (lres[-2][0], fres[3][0], '1 ;cycle')  # C-HCA
    pair7 = '%5s %5s %12s' % (lres[-2][0], fres[4][0], '1 ;cycle')  # C-_CA
    pair8 = '%5s %5s %12s' % (lres[-2][0], fres[-2][0], '1 ;cycle')  # C-C (carbonyls)
    pair9 = '%5s %5s %12s' % (lres[-1][0], fres[1][0], '1 ;cycle')  # O-HN
    pair10 = '%5s %5s %12s' % (lres[-1][0], fres[2][0], '1 ;cycle')  # O-CA


    # angles
    angle1 = '%5s %5s %5s %12s' % (lres[2][0], lres[-2][0], fres[0][0], '5 ;cycle')
    angle2 = '%5s %5s %5s %12s' % (lres[-1][0], lres[-2][0], fres[0][0], '5 ;cycle')
    angle3 = '%5s %5s %5s %12s' % (lres[-2][0], fres[0][0], fres[1][0], '5 ;cycle')
    angle4 = '%5s %5s %5s %12s' % (lres[-2][0], fres[0][0], fres[2][0], '5 ;cycle')

    # dihedrals
    died1 = '%5s %5s %5s %5s %12s' % (lres[0][0], lres[2][0], lres[-2][0], fres[0][0], '9 ;cycle')
    died2 = '%5s %5s %5s %5s %12s' % (lres[3][0], lres[2][0], lres[-2][0], fres[0][0], '9 ;cycle')
    died3 = '%5s %5s %5s %5s %12s' % (lres[4][0], lres[2][0], lres[-2][0], fres[0][0], '9 ;cycle')
    died4 = '%5s %5s %5s %5s %12s' % (lres[2][0], lres[-2][0], fres[0][0], fres[1][0], '9 ;cycle')
    died5 = '%5s %5s %5s %5s %12s' % (lres[2][0], lres[-2][0], fres[0][0], fres[2][0], '9 ;cycle')
    died6 = '%5s %5s %5s %5s %12s' % (lres[-1][0], lres[-2][0], fres[0][0], fres[1][0], '9 ;cycle')
    died7 = '%5s %5s %5s %5s %12s' % (lres[-1][0], lres[-2][0], fres[0][0], fres[2][0], '9 ;cycle')
    died8 = '%5s %5s %5s %5s %12s' % (lres[-2][0], fres[0][0], fres[2][0], fres[3][0], '9 ;cycle')
    died9 = '%5s %5s %5s %5s %12s' % (lres[-2][0], fres[0][0], fres[2][0], fres[4][0], '9 ;cycle')
    died10 = '%5s %5s %5s %5s %12s' % (lres[-2][0], fres[0][0], fres[2][0], fres[-2][0], '9 ;cycle')

    # dihedrals
    impr1 = '%5s %5s %5s %5s %12s' % (lres[-2][0], lres[2][0], fres[0][0], lres[-1][0], '2 ;cycle')
    impr2 = '%5s %5s %5s %5s %12s' % (fres[0][0], lres[-2][0], fres[2][0], fres[1][0], '2 ;cycle')

    # write the new topology file
    g = open(args.output, 'w')

    # Write prologue
    for line in paragraphs['prologue']:
        g.write(line)
    del par_order_in[0]
    del paragraphs['prologue']

    # Write other blocks like `[ system ]` `[ moleculetype ]` `#include`
    for parname in par_order_in:
        if parname.startswith('#include') or parname.startswith('['):
            if '[ molecules ]' in parname:  # Gromacs wants this at the end so it can easilty append number of molecules and ions added using `solvate` and `genion`
                epilogue = paragraphs[parname]
            else:
                for line in paragraphs[parname]:
                    g.write(line)
                g.write('\n')
            del paragraphs[parname]


    for i, line in enumerate(l):
        d = l[i].split()
        if i < prologue_end:
            continue
        if not l[i].split():
            continue
        if i == b[0]:
            g.write('\n')
        if d[0] == ';' or d[0] == '[':
            g.write(l[i])
        elif d[0].startswith('Protein'):
            g.write(l[i]+'\n')
        elif d[0].isdigit() and d[1].isupper():
            g.write(l[i])
        elif d[0].isdigit() and d[1].isdigit():
            g.write(l[i])
        if d[0] == lres[-2][0] and d[1] == lres[-1][0] and d[2] == '1':  # Add the new bond as the last one
            g.write(bond + '\n')
            g.write(' ' + '\n')
        elif d[0] == lres[0][0] and d[1] == lres[-1][0] and d[2] == '1':
            g.write(pair1 + '\n')
        elif d[0] == lres[2][0] and d[1] == lres[-3][0] and d[2] == '1':
            g.write(pair2 + '\n')
            g.write(pair3 + '\n')
        elif d[0] == lres[3][0] and d[1] == lres[-1][0] and d[2] == '1':
            g.write(pair4 + '\n')
        elif d[0] == lres[4][0] and d[1] == lres[-1][0] and d[2] == '1':
            g.write(pair5 + '\n')
        if d == l[int(q[0])-2].split():
            g.write(pair6 + '\n')
            g.write(pair7 + '\n')
            g.write(pair8 + '\n')
            g.write(pair9 + '\n')
            g.write(pair10 + '\n')
            g.write(' ' + '\n')
        if d[0] == lres[2][0] and d[1] == lres[-2][0] and d[2] == lres[-1][0] and d[3] == '5':
            g.write(angle1 + '\n')
            g.write(angle2 + '\n')
            g.write(angle3 + '\n')
            g.write(angle4 + '\n')
            g.write(' ' + '\n')
        if d[0] == lres[0][0] and d[1] == lres[2][0] and d[2] == lres[-2][0] and d[3] == lres[-1][0] and d[4] == '9':
            g.write(died1 + '\n')
        elif d[0] == lres[3][0] and d[1] == lres[2][0] and d[2] == lres[-2][0] and d[3] == lres[-1][0] and d[4] == '9':
            g.write(died2 + '\n')
        elif d[0] == lres[4][0] and d[1] == lres[2][0] and d[2] == lres[-2][0] and d[3] == lres[-1][0] and d[4] == '9':
            g.write(died3 + '\n')
        if d == l[int(w[1])-2].split():
            g.write(died4 + '\n')
            g.write(died5 + '\n')
            g.write(died6 + '\n')
            g.write(died7 + '\n')
            g.write(died8 + '\n')
            g.write(died9 + '\n')
            g.write(died10 + '\n')
            g.write(' ' + '\n')
        if d == l[int(z[0])-2].split():
            g.write(impr1 + '\n')
            g.write(impr2 + '\n')
            g.write(' ' + '\n')
        if 'cmap' in d:
            g.write(';  ai    aj    ak    al    am funct'+'\n')
            for i in range(0, len(C)):
                map = ('%5s %5s %5s %5s %5s %5s \n' %
                       (C[i-1], N[i], CA[i], C[i], N[i-(len(C)-1)], '1'))
                g.write(map)
            g.write('\n')
            for i, par in enumerate(paragraphs):
                for line in paragraphs[par]:
                    g.write(line)
                g.write('\n')
            for line in epilogue:  # Writes the `[ molecules ]` block, which gromacs wants to be at the end of the file
                g.write(line)
            break
    print(f"New topology was stored in `{args.output}`.")
    g.close()

if __name__ = '__main__':
    main()
