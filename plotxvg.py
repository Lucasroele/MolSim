import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import numpy as np
import itertools


def getFileNames(extension):
    # List all files in the current working directory
    current_directory = os.getcwd()
    files_in_directory = os.listdir(current_directory)
    # Filter out directories from the list
    filenames = [file for file in files_in_directory \
                 if os.path.isfile(os.path.join(current_directory, file)) \
                 and file.endswith(extension)]
    filenames.sort()
    return filenames


def getDataFromXmgrFiles(_filenames):
    """
    Assumes that strings are double quote enclosed
    Returns
    dict with
        <filename> np.arrays with columns in xmgrace stored as rows
        <cols> max#ofcols
        <y_axis_labels> nested list
    """
    if isinstance(_filenames, str):
        filenames = [_filenames]
    else:
        filenames = _filenames
    # Initialize lists to store the data
    data = {filenames[i]:[] for i in range(len(filenames))}
    y_axis_labels = {}
    cols = 0  # Number of columns used for plotting
    # Read data from each file
    for i, filename in enumerate(filenames):
        y_axis_labels[filename] = []
        with open(filename, 'r') as file:
            j = 0
            for line in file:
                if line.startswith('#'):
                    continue
                if line.strip() == '':
                    continue
                if line.startswith('@'):  # parse in column names
                    infoline = line.split()
                    if len(infoline) >= 4 and infoline[2] == "legend" and infoline[1][0] == "s" \
                    or len(infoline) >= 5 and infoline[1] == "legend" and infoline[2] == "string":
                        if infoline[-1].endswith('"') and not infoline[-1].startswith('"'):
                            label = infoline[-1]
                            for substr in infoline[-2:0:-1]:
                                label = substr + " " + label
                                if substr.startswith('"'):
                                    break
                        else:
                            label = infoline[-1]
                        label = label.strip('"')
                        y_axis_labels[filename].append(label)
                    continue
                infoline = line.split()
                data[filename].append([])
                for k, col in enumerate(infoline):
                    data[filename][j].append(float(col))
                if k > cols:
                    cols = k
                j += 1

    # When no legend in the xvg present
    for i, filename in enumerate(filenames):
        if len(y_axis_labels[filename]) != len(data[filename][0]) - 1:
            y_axis_labels[filename] = list(range(len(data[filename])+1)[1::])

    # Create the subplot figure
    data = {filenames[i]:np.array(data[filenames[i]]).T for i in range(len(data))}
    data["cols"] = cols
    data["y_axis_labels"] = y_axis_labels
    return data


def makeFig(filenames, data):
    # Create 2d list with plot titles
    subplot_titles = [[filename if i < len(data[filename]) - 1 else "" for i in range(data["cols"])] for filename in filenames]
    fig = make_subplots(rows=len(filenames), cols=data["cols"], subplot_titles = list(itertools.chain(*subplot_titles)))
    # Add subplots
    for i, filename in enumerate(filenames):
        print(filename)
        print(data["y_axis_labels"][filename])
        for j, y in enumerate(data[filename]):
            if j == 0:
                continue
            fig.add_trace(go.Scatter(x = data[filename][0],
                                     y = y),
                                     row = i + 1,
                                     col = j)
            fig.update_xaxes(title_text="Time",
                             row = i + 1,
                             col = j)
            fig.update_yaxes(title_text = data["y_axis_labels"][filename][j - 1],
                             row = i + 1,
                             col = j)

    # Update layout
    fig.update_layout(height=300 * len(filenames), showlegend=False)
    return fig


def main():
    filenames = getFileNames(".xvg")
    data = getDataFromXmgrFiles(filenames)
    fig = makeFig(filenames, data)
    # Display the plot
    fig.show()


main()
