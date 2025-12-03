"""
Git test test2
"""

import argparse
import itertools
import os
import numpy as np
from pathlib import Path

from molsim.parsers import XvgParser
#from molsim.utils import getFileNames


def register(subparsers):
    parser = subparsers.add_parser('frame_times',
                                    help='Finds the frame times of frames with an observed set of r values that can be specified in several ways. The program assumes the second column of the .xvg file contains ascending r values over time.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)



def parseArguments():
    parser = argparse.ArgumentParser(prog='frame_times.py',
                                     description='Finds the frame times of frames with an observed set of r values that can be specified in several ways. The program assumes the second column of the .xvg file contains ascending r values over time.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('xvgfile',
                        type=str,
                        help='filename of the .xvg with r values')           # positional argument
    parser.add_argument('-b',
                        '--begin',
                        type=float,
                        help='value of the first frame')           # positional argument
    parser.add_argument('-e',
                        '--end',
                        type=float,
                        help='value of the last frame')           # positional argument
    parser.add_argument('-nf',
                        '--n_frames',
                        type=int,
                        help='number of frames')           # positional argument
    parser.add_argument('-dr',
                        '--dr',
                        type=float,
                        help='spacing between separate r values')
    parser.add_argument('-rarr',
                        '--rarray',
                        nargs='+',
                        type=float,
                        help='space separated r values')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='write the output to a file instead of the terminal.')
    parser.add_argument('-log',
                         '--log',
                         type=str,
                         default='distances.dat',
                         help='specify a file different from `distances.dat` to output logging into')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='overwrite any files already present')
    return parser


def validateArguments(args):
    # It doesn't make sense to specify some combinations of arguments
    if args.rarray is not None:
        for i in ['dr', 'n_frames', 'end', 'begin']:
            assert getattr(args, i) is None, f"It doesn't make sense to specify both `rarray` and `{i}`."
    else:
        assert args.begin is not None, "The program doesn't know how to determine the first frame. If `rarray` is not specified, `begin` must be specified."
        assert args.end is not None or args.n_frames is not None or args.dr is not None, "The program doesn't know how to determine the r values."
        if args.dr is not None:
            assert args.n_frames is not None or args.end is not None, "If `dr` is specified, either `n_frames` or `end` must be specified."
            assert not (args.n_frames is not None and args.end is not None), "It doesn't make sense to specify `begin`, `dr`, `n_frames` and `end` together."
        if args.end is not None:
            assert args.end > args.begin, "`end` must be larger than `begin`."
    assert Path(args.xvgfile).is_file(), f"`{args.xvgfile}` does not exist."
    if not args.force:
        assert args.output is None or not os.path.isfile(args.output), f"`{args.output}` already exists, specify `force` to ignore this check and overwrite it"
        assert not Path(args.log).is_file(), f"{args.log} already exists, specify `force` to ignore this check and overwrite it"
    return args


def main(args):
    args = validateArguments(args)
    file_index = 0
    xvgObj = XvgParser(args.xvgfile)
    n_frames_in_file = len(xvgObj.data[file_index][0])
    minmax = []
    for i, row in enumerate(xvgObj.data[file_index]):
        minmax.append([])
        minmax[i].append(min(row))
        minmax[i].append(max(row))
    total_dist = minmax[1][1] - minmax[1][0]
    badargs_mes = "Error: The arguments do not fit the input file."
    rarray = []  # Becomes a np array
    dr = None
    # Filling rarray if not specified
    if args.rarray is None:
        assert args.begin > minmax[1][0] and args.begin < minmax[1][1], badargs_mes
        rarray.append(args.begin)
        if args.end is not None:
            assert args.end < minmax[1][1], badargs_mes
            if args.n_frames is not None:
                assert args.n_frames < n_frames_in_file, badargs_mes
                dr = (args.end - args.begin) / float(args.n_frames-1)
                for i in range(1, args.n_frames - 1):
                    rarray.append(args.begin + i*dr)
            elif args.dr is not None:
                assert int((args.end - args.begin) / args.dr) + 1< n_frames_in_file, badargs_mes
                # /\ could do better maths here /\
                for i in itertools.count():
                    assert args.begin < args.end and args.dr > 0, "Infinite loop detected, check your arguments."
                    rarray.append(args.begin + i*args.dr)
                    if rarray[-1] + args.dr > args.end:
                        break
            rarray.append(args.end)
        elif args.dr is not None and args.n_frames is not None:
            assert args.dr * args.n_frames < total_dist, badargs_mes
            for i in range(1, args.n_frames):
                rarray.append(args.begin + i*args.dr)
        else:
            assert False, "Something went wrong."
        rarray = np.array(rarray)
    else:
        rarray = np.array(args.rarray)
    ###### Print what you will be doing into log file
    output = []
    for i, curval in enumerate(rarray):
        tempout = 0
        timeout = None
        for j, val in enumerate(xvgObj.data[file_index][1][-1:0:-1]):
            if abs(curval - tempout) > abs(curval - val):
                tempout = val
                timeout = xvgObj.data[file_index][0][-1-j]
        output.append(str(int(timeout)))
    print(' '.join(output))


if __name__ == "__main__":
    args = parseArguments()
    main(args)
