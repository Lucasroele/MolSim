"""
Git test test2
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots



class XvgParser:
    """
    Usage:
    initialized using a list of (.xvg) filenames or a single filename

    Holds:
    filenames = [filename_s]
    data = [np.ndarray]         # xvg-cols in rows
    metadata = [{key, val}]
    size = number of entries
    """

    def __init__(self, filename_s, path=None):
        self.data = []
        self.metadata = []
        if isinstance(filename_s, str):
            self.filenames = [filename_s]
        elif isinstance(filename_s, list):
            self.filenames = filename_s
        else:
            raise ValueError("Input must be a string or a list of strings")
        if path is not None:
            assert os.path.exists(path)
            assert os.path.isdir(path)
            if not path.endswith('/'):
                path += '/'
        else:
            path = ''
        for i, filename in enumerate(self.filenames):
            assert filename.endswith('.xvg')
            assert os.path.exists(path + filename)
            assert os.path.isfile(path + filename)
            self.filenames[i] = path + filename
        for i, filename in enumerate(self.filenames):
            reading_data = False
            tempdata = []
            self.metadata.append({})
            self.metadata[i]['comments'] = []
            self.metadata[i]['columns'] = []
            with open(filename, 'r') as file:
                prev_line = ''
                for line in file:
                    if line.startswith('#'):
                        self.metadata[i]['comments'].append(line)
                        prev_line = line
                        continue
                    if len(line.strip(' ')) == 0:
                        prev_line = line
                        continue
                    if line.startswith('@'):  # parse in column names
                        infoline = line.split()
                        if len(infoline) <= 2:
                            prev_line = line
                            continue
                        # infolinen will have at least 3 elements
                        if infoline[1] == 'title':
                            self.metadata[i]['title'] = self.findString(
                                infoline[2:])
                        elif infoline[1] == 'xaxis':
                            if infoline[2] == 'label' and len(infoline) > 3:
                                self.metadata[i]['xaxis'] = self.findString(
                                    infoline[3:])
                            else:
                                self.metadata[i]['xaxis'] = self.findString(
                                    infoline[2:])
                        elif infoline[1] == 'yaxis':
                            if infoline[2] == 'label' and len(infoline) > 3:
                                self.metadata[i]['yaxis'] = self.findString(
                                    infoline[3:])
                            else:
                                self.metadata[i]['yaxis'] = self.findString(
                                    infoline[2:])
                        elif len(infoline) >= 4 and infoline[2] == "legend" and infoline[1][0] == "s":
                            self.metadata[i]['columns'].append(
                                self.findString(infoline[3:]))
                        elif len(infoline) >= 5 and infoline[1] == "legend" and infoline[2] == "string":
                            self.metadata[i]['columns'].append(
                                self.findString(infoline[4:]))
                        continue
                    # Reading data lines now
                    infoline = line.split()
                    if not reading_data:
                        reading_data = True
                        for cel in infoline:
                            tempdata.append([])
                        if len(self.metadata[i]['columns']) != len(tempdata):
                            if prev_line == self.metadata[i]['comments'][-1] and prev_line.split() == len(tempdata):
                                self.metadata[i]['columns'] = [
                                    cel for cel in prev_line.split()]
                            elif len(self.metadata[i]['columns']) == 1:
                                self.metadata[i]['columns'].insert(0, '#')
                            else:
                                self.metadata[i]['columns'] = [
                                    str(leg) for leg, _ in enumerate(tempdata)]
                    for k, cel in enumerate(infoline):
                        tempdata[k].append(float(cel))
                    prev_line = line
                self.data.append(np.array(tempdata))
        self.size = len(self.data)

    def findString(self, list_of_strings):
        """
        ['"a', 'b', 'c"'] -> "a b c"
        """
        if len(list_of_strings) == 1:
            return list_of_strings[0].strip('"')
        ret_string = ''
        found_start = False
        for string in list_of_strings:
            if not found_start and string.startswith('"'):
                found_start = True
                ret_string += string
                continue
            elif found_start and string.endswith('"'):
                ret_string += ' ' + string
                break
            elif found_start:
                ret_string += ' ' + string
        assert len(ret_string) != 0
        return ret_string.strip('"')

    def indexOfFilename(self, filename):
        for i, filename_ in enumerate(self.filenames):
            if filename == filename_:
                return i

    def indexIsRetrievable(self, index):
        if isinstance(index, int):
            assert index in range(self.size)
        elif isinstance(index, str):
            assert index in self.filenames
        else:
            raise ValueError(f"{index} is not retrievable")
        return True

    def get_numpy(self, index):  # returns as rows of datapoints
        assert self.indexIsRetrievable(index)
        if isinstance(index, str):
            return self.data[self.indexOfFilename(index)]
        else:  # Must be int
            return self.data[index]

    def get_pandas(self, index):  # returns as columns of datapoints
        assert self.indexIsRetrievable(index)
        if isinstance(index, str):
            print(self.metadata[self.indexOfFilename(index)]['columns'])
            return pd.DataFrame(data=self.data[self.indexOfFilename(index)].T,
                                columns=self.metadata[self.indexOfFilename(index)]['columns'])
        else:  # Must be int
            return pd.DataFrame(data=self.data[index],
                                columns=self.metadata[index]['columns'])

    def plot(self, all_in_one=False, subplot_titles=None, showlegend=True):
        """
        all_in_one:     overlay all data into a single plot
        subplot_titles: define titles for the plots (default = filenames)
        """
        colors = px.colors.qualitative.Alphabet
        if all_in_one:
            rows = 1
            if subplot_titles is not None:
                subplot_titles = [str(subplot_titles)]
            else:
                subplot_titles = ["All together"]
        else:
            rows = self.size
            if subplot_titles is not None:
                assert len(subplot_titles) == len(self.filenames), f"The argument `subplot_titles` should be of lenght: {self.size}"
            else:
                subplot_titles = self.filenames
        fig = make_subplots(rows=rows,
                            cols=1,
                            subplot_titles=subplot_titles)
        for i, data in enumerate(self.data):  # subplots
            for j, row in enumerate(data):  # rows
                if j == 0:
                    continue
                if all_in_one:
                    fig.add_trace(go.Scatter(x=data[0],
                                             y=row,
                                             legendgroup=self.filenames[i] + self.metadata[i]['columns'][j],
                                             name=self.metadata[i]['columns'][j]),
                                  row=1,
                                  col=1)
                else:
                    use_row = i + 1
                    trace = go.Scatter(x=data[0],
                                       y=row,
                                       line=dict(color=colors[j]),
                                       legendgroup=self.metadata[i]['columns'][j],
                                       name=self.metadata[i]['columns'][j])
                    if i > 0:
                        trace['showlegend'] = False
                    fig.add_trace(trace,
                                  row=use_row,
                                  col=1)
            if rows > 1 or all_in_one and i == 0:
                if 'xaxis' in self.metadata[i]:
                    fig.update_xaxes(title_text=self.metadata[i]['xaxis'],
                                     row=i + 1,
                                     col=1)
                if 'yaxis' in self.metadata[i]:
                    fig.update_yaxes(title_text=self.metadata[i]['yaxis'],
                                     row=i + 1,
                                     col=1)
        # Update layout
        fig.update_layout(height=300 * len(self.filenames), showlegend=showlegend)
        fig.show()


def parseArguments():
    parser = argparse.ArgumentParser(prog='frameTimesGetter.py',
                                     description='Finds the frame times of frames with an observed set of r values that can be specified in several ways.',
                                     epilog='Written by Lucas Roeleveld')

    parser.add_argument('filename1',
                        nargs=1,
                        help='filename of the .xvg with r values')           # positional argument
    parser.add_argument('-b',
                        '--begin',
                        nargs=1,
                        type=float,
                        help='value of the first frame')           # positional argument
    parser.add_argument('-e',
                        '--end',
                        nargs=1,
                        type=float,
                        help='value of the last frame')           # positional argument
    parser.add_argument('-nf',
                        '--n_frames',
                        nargs=1,
                        type=int,
                        help='number of frames')           # positional argument
    parser.add_argument('-dr',
                        '--dr',
                        nargs=1,
                        type=float,
                        help='spacing between separate r values')
    parser.add_argument('-rarr',
                        '--rarray',
                        nargs='+',
                        type=float,
                        help='space separated r values')
    parser.add_argument('-o',
                        '--output',
                        nargs=1,
                        type=str,
                        help='write the output to a file instead of the terminal.')
    parser.add_argument('-log',
                         '--log',
                         nargs=1,
                         type=str,
                         default='distances.dat',
                         help='specify a file different from `distances.dat` to output logging into')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='overwrite any files already present')
    args = parser.parse_args()
    # It doesn't make sense to specify some combinations of arguments
    if args.rarray is not None:
        for i in ['dr', 'n_frames', 'end', 'begin']:
            assert getattr(args, i) is None, f"It doesn't make sense to specify both --rarray/-rarr and {i}."
    if args.dr is not None:
        assert args.n_frames is None, "It doesn't make sense to specify both --dr/-dr and --n_frames/-nf."
    assert os.path.isfile(args.filename1), f"{args.filename1} does not exist"
    if not args.force:
        assert args.output is None  or not os.path.isfile(args.output), f"{args.output} already exists, specify `-f` to ignore this check and overwrite it"
        assert not os.path.isfile(args.log), f"{args.log} already exists, specify `-f` to ignore this check and overwrite it"
    return args


def getFileNames(extension):
    # List all files in the current working directory
    current_directory = os.getcwd()
    files_in_directory = os.listdir(current_directory)
    # Filter out directories from the list
    filenames = [file for file in files_in_directory
                 if os.path.isfile(os.path.join(current_directory, file))
                 and file.endswith(extension)]
    filenames.sort()
    return filenames


def main():
    badargs_mes = "You should specify either an array of values or use the other options to specify it indirectly."
    MAX_FRAMES = 1000
    args = parseArguments()
    if args.rarray is None:
        rarray = []  # Becomes a np array
        if args.begin is None:  # We don't do no backwards specifiying from args.end
            sys.exit(badargs_mes)
        else:
            rarray.append(args.begin)
            if args.end is not None and args.n_frames is not None:
                assert args.n_frames < MAX_FRAMES
                assert args.end > args.begin
                dr = (args.end - args.begin) / float(n_frames-1)
                for i in range(args.n_frames-1)[1:]:
                    rarray.append(args.begin + i*dr)
            elif args.end is not None and args.dr is not None:
                assert args.end > args.begin
                for i in range(MAXFRAMES)[1:]:
                    if args.begin + i*args.dr < args.end:
                        rarray.append(args.begin + i*args.dr)
                    else:
                        break
            elif args.n_frames is not None and args.dr is not None:
                assert args.n_frames < MAX_FRAMES
                for i in range(args.n_frames)[1:]:
                    rarray.append(args.begin + i*args.dr)
            else:
                sys.exit(badargs_mes)
        rarray = np.array(rarray)
    else:
        rarray = np.array(args.rarray)
    startval = np.min(rarray)
    endval = np.max(rarray)

    xvgObj = XvgParser(args.filename1)
    index = 0

    minmax = []
    for i, row in enumerate(xvgObj.data[index]):
        minmax.append([])
        minmax[i].append(min(row))
        minmax[i].append(max(row))

    assert args.startval > minmax[1][0] and endval < minmax[1][1], f"The specified range is not covered by the values inside {args.filename1}."

    ###### Print what you will be doing into log file
    ### use rarray instead of equidistant spacing

    output = []
    for i in range(args.n_frames):
        tempout = 0
        timeout = None
        curval = args.startval + i*args.stepval
        for j, val in enumerate(xvgObj.data[index][1][-1:0:-1]):
            if abs(curval - tempout) > abs(curval - val):
                tempout = val
                timeout = xvgObj.data[index][0][-1-j]
        output.append(str(int(timeout)))
    print(' '.join(output))
main()
