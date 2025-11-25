import sys
import os
import numpy as np

from ..parsers.xvgParser import XvgParser
from ..utils import getFileNames


# CLI registration as subcommand
def register(subparsers):
    parser = subparsers.add_parser('xvg_min',
                                   help='Print the smallest and largest numbers in each column of an .xvg file')
    parser.add_argument('filename1', # positional argument
                        nargs='?',
                        default=None,
                        help='the .xvg file.')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='append the output to a file instead of the terminal.')
    parser.set_defaults(func=main)


def getDataFromXmgrFiles(_filenames):
    """
    Assumes that strings are double quote enclosed
    Returns
    dict with
        <filename> np.arrays with columns in xmgrace stored as rows
        <cols> max#ofcols
        <legend_entries> nested list
    """
    if isinstance(_filenames, str):
        filenames = [_filenames]
    else:
        filenames = _filenames
    # Initialize lists to store the data
    data = {filenames[i]: [] for i in range(len(filenames))}
    legend_entries = {}
    cols = 0  # Number of columns used for plotting
    # Read data from each file
    reading_cols = False
    prev_line = ''
    for i, filename in enumerate(filenames):
        legend_entries[filename] = []
        with open(filename, 'r') as file:
            j = 0
            for m, line in enumerate(file):
                if line.startswith('#'):
                    prev_line = line
                    continue
                if line.strip() == '':
                    prev_line = line
                    continue
                if line.startswith('@'):  # parse in column names
                    infoline = line.split()
                    if len(infoline) >= 4 and infoline[2] == "legend" and infoline[1][0] == "s" \
                    or len(infoline) >= 5 and infoline[1] == "legend" and infoline[2] == "string":
                        label = ""
                        if infoline[-1].endswith('"') and not infoline[-1].startswith('"'):
                            for substr in reversed(infoline):
                                label = substr + " " + label
                                if substr.startswith('"'):
                                    break
                        else:
                            label = infoline[-1]
                        label = label.strip('"')
                        legend_entries[filename].append(label)
                    prev_line = line
                    continue
                # Reading data lines now
                infoline = line.split()
                assert len(infoline) > 1, "This script assumes the x axis is defined in the .xvg file"
                data[filename].append([])
                for k, col in enumerate(infoline):
                    data[filename][j].append(float(col))
                if k > cols:
                    cols = k
                j += 1
                # Make sure the legend_entries are good
                if not reading_cols:
                    if len(legend_entries[filename]) < len(infoline) - 1:
                        if len(prev_line.split()) == len(infoline):
                            legend_entries[filename] = [prev_line.split()[1]]
                            for name in prev_line.split()[2:]:
                                legend_entries[filename].append(name)
                        else:
                            legend_entries[filename] = [str(a) for a in range(len(infoline) + 1)[1:]]


                prev_line = line
                reading_cols = True
            assert len(legend_entries[filename]) == len(data[filename][0]) - 1, "legend entries not parsed well"

    # Create the subplot figure
    data = {filenames[i]: np.array(data[filenames[i]]).T for i in range(len(data))}
    data["cols"] = cols
    data["legend_entries"] = legend_entries
    return data


def main(args):
    #args = parseArguments()
    if args.filename1 is None:
        args.filename1 = getFileNames(".xvg")[0]
    running_folder = os.getcwd()
    if not os.path.exists(args.filename1):
        sys.exit("Error: " + "The file `" + running_folder + "/" + args.filename1 + "` does not exist.")

    #data = getDataFromXmgrFiles(args.filename1)
    xvgObj = XvgParser(args.filename1)
    index = 0

    minmax = []
    for i, row in enumerate(xvgObj.data[index]):
        minmax.append([])
        minmax[i].append(min(row))
        minmax[i].append(max(row))

    headers = [] 
    if len(minmax) == 2 and 'xaxis' in xvgObj.metadata[index] and 'yaxis' in xvgObj.metadata[index]:
        headers = [xvgObj.metadata[index]['xaxis'], xvgObj.metadata[index]['yaxis']]

    if args.output is not None:
        if os.path.exists(args.output):
            writemode = 'a'
        else:
            writemode = 'w'
        with open(args.output, writemode) as file:
            file.write(f"\nFetched from {args.filename1}:\n")
            for i, row in enumerate(xvgObj.data[index]):
                if len(headers) > i:
                    file.write(f"{headers[i]}:\n")
                else:
                    file.write(f"{xvgObj.metadata[index]['columns'][i]}:\n")
                file.write(f"\tmin: {minmax[i][0]:.3e}\n")
                file.write(f"\tmax: {minmax[i][1]:.3e}\n")

    else:
        print(f"\nFetched from {args.filename1}:")
        for i, row in enumerate(xvgObj.data[index]):
            if len(headers) > i:
                print(f"{headers[i]}:\n")
            else:
                print(f"{xvgObj.metadata[index]['columns'][i]}:\n")
            print(f"\tmin: {minmax[i][0]:.3e}\n")
            print(f"\tmax: {minmax[i][1]:.3e}\n")
            #print(f"{xvgObj.metadata[index]['columns'][1]}:\n\tmin: {minmax[i - 1][0]:.3e}\n\tmax: {minmax[i - 1][1]:.3e}")

if __name__ == "__main__":
    main()
