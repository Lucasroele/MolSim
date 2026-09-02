import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import itertools
import argparse
from pathlib import Path

from molsim.utils.files import getFileNames
from molsim.parsers.xvgParser import XvgParser


def register(subparsers):
    parser = subparsers.add_parser('plot_xvg',
                                   help='Plot an .xvg file or all .xvg files in a directory using plotly. If no argument is supplied the xvg files in the working directory will be plotted.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='plot_xvg',
                                     description='Plot an .xvg file or all .xvg files in a directory using plotly. If no argument is supplied the xvg files in the working directory will be plotted.',
                                     epilog='Written by Lucas Roeleveld')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('filenames',
                        nargs='*',
                        help='the .xvg file(s).')           # positional argument
    parser.add_argument('-d',
                        '--dir',
                        type=str,
                        default=None,
                        help='the directory with .xvg files to be plotted. Only makes sense if no filenames are supplied.')  # optional argument
    return parser


def validateArguments(args):
    running_folder = os.getcwd()
    if args.dir is not None and args.filenames is not []:
        print("Note: The --dir argument is ignored when filenames are supplied.")
    if args.filenames == []:
        if args.dir is None:
            args.filenames = getFileNames(".xvg")
        else:
            assert Path(args.dir).is_dir(), "Error: The directory `" + args.dir + "` does not exist."
            args.filenames = getFileNames(".xvg", path=args.dir)
        assert len(args.filenames) > 0, "Error: " + "No .xvg files found in " + (args.dir if args.dir is not None else running_folder) + "."
    else:
        for filename in args.filenames:
            assert os.path.exists(filename), "Error: The file `" + filename + "` does not exist."
    return args


def makeFig(filenames, dataset):
    # Create 2d list with plot titles
    cols = max([len(dataset[filename]) - 1 for filename in filenames])
    subplot_titles = [["" if j > len(dataset[filename]) - 1 else dataset.metadata[i]['columns'][j] for j in range(1+cols)[1:]] for i, filename in enumerate(filenames)]
    fig = make_subplots(rows=len(filenames), cols=cols, subplot_titles = list(itertools.chain(*subplot_titles)))
    # Add subplots
    for i, filename in enumerate(filenames):
        for j, y in enumerate(dataset[filename]):
            if j == 0:
                continue
            fig.add_trace(go.Scatter(x = dataset[filename][0],
                                     y = y),
                                     row = i + 1,
                                     col = j)
            fig.update_xaxes(title_text=dataset.metadata[i]['xaxis'],
                             row = i + 1,
                             col = j)
            fig.update_yaxes(title_text = dataset.metadata[i]['yaxis'],
                             row = i + 1,
                             col = j)

    # Update layout
    fig.update_layout(height=300 * len(filenames), showlegend=False)
    return fig


def main(args):
    args = validateArguments(args)
    if args.dir is not None:
        dataset = XvgParser(args.filenames, path=args.dir)
    else:
        dataset = XvgParser(args.filenames)
    #filenames = getFileNames(".xvg")
    #dataset = getDataFromXmgrFiles(args.filenames)
    fig = makeFig(args.filenames, dataset)
    # Display the plot
    fig.show()



if __name__ == "__main__":
    args = parseArguments()
    main(args)
