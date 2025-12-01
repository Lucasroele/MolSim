import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import itertools
import argparse
from pathlib import Path

from molsim.utils import getFileNames
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
                        nargs='?',
                        default=None,
                        help='the .xvg file(s).')           # positional argument
    parser.add_argument('-d',
                        '--dir',
                        type=str,
                        default=None,
                        help='the directory with .xvg files to be plotted. Only makes sense if no filenames are supplied.')  # optional argument
    return parser


def validateArguments(args):
    running_folder = os.getcwd()
    if args.dir is not None and args.filenames is not None:
        print("Note: The --dir argument is ignored when filenames are supplied.")
    if args.filenames is None:
        if args.dir is None:
            args.filenames = getFileNames(".xvg")
        else:
            assert Path(args.dir).is_dir(), "Error: The directory `" + args.dir + "` does not exist."
            args.filenames = getFileNames(".xvg", path=args.dir)
        assert len(args.filenames) > 0, "Error: " + "No .xvg files found in " + (args.dir if args.dir is not None else running_folder) + "."
    elif isinstance(args.filenames, str):
        assert os.path.exists(args.filename1), "Error: The file `" + args.filenames + "` does not exist."
    elif isinstance(args.filenames, list):
        for filename in args.filenames:
            assert os.path.exists(filename), "Error: The file `" + filename + "` does not exist."
    return args

#def getDataFromXmgrFiles(_filenames):
#    """
#    Assumes that strings are double quote enclosed
#    Returns
#    dict with
#        <filename> np.arrays with columns in xmgrace stored as rows
#        <cols> max#ofcols
#        <y_axis_labels> nested list
#    """
#    if isinstance(_filenames, str):
#        filenames = [_filenames]
#    else:
#        filenames = _filenames
#    # Initialize lists to store the data
#    data = {filenames[i]:[] for i in range(len(filenames))}
#    y_axis_labels = {}
#    cols = 0  # Number of columns used for plotting
#    # Read data from each file
#    for i, filename in enumerate(filenames):
#        y_axis_labels[filename] = []
#        with open(filename, 'r') as file:
#            j = 0
#            for line in file:
#                if line.startswith('#'):
#                    continue
#                if line.strip() == '':
#                    continue
#                if line.startswith('@'):  # parse in column names
#                    infoline = line.split()
#                    if len(infoline) >= 4 and infoline[2] == "legend" and infoline[1][0] == "s" \
#                    or len(infoline) >= 5 and infoline[1] == "legend" and infoline[2] == "string":
#                        if infoline[-1].endswith('"') and not infoline[-1].startswith('"'):
#                            label = infoline[-1]
#                            for substr in infoline[-2:0:-1]:
#                                label = substr + " " + label
#                                if substr.startswith('"'):
#                                    break
#                        else:
#                            label = infoline[-1]
#                        label = label.strip('"')
#                        y_axis_labels[filename].append(label)
#                    continue
#                infoline = line.split()
#                data[filename].append([])
#                for k, col in enumerate(infoline):
#                    data[filename][j].append(float(col))
#                if k > cols:
#                    cols = k
#                j += 1
#
#    # When no legend in the xvg present
#    for i, filename in enumerate(filenames):
#        if len(y_axis_labels[filename]) != len(data[filename][0]) - 1:
#            y_axis_labels[filename] = list(range(len(data[filename])+1)[1::])
#
#    # Create the subplot figure
#    data = {filenames[i]:np.array(data[filenames[i]]).T for i in range(len(data))}
#    data["cols"] = cols
#    data["y_axis_labels"] = y_axis_labels
#    return data



def makeFig(filenames, dataset):
    # Create 2d list with plot titles
    cols = max([len(dataset[filename]) - 1 for filename in filenames])
    subplot_titles = [[filename if i < len(dataset[filename]) - 1 else "" for i in range(cols)] for filename in filenames]
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
