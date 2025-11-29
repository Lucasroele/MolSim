import os
import numpy as np
import pandas as pd
import polars
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


class XvgParser:
    """
    Instantiation:
        xvg = XvgParser(filename_s, path=None)
    
    filenames = [filename_s]
    data = [np.ndarray]         # Cols in rows
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
            assert os.path.exists(path), f"Path `{path}` does not exist"
            assert os.path.isdir(path)
            if not path.endswith('/'):
                path += '/'
        else:
            path = ''
        for i, filename in enumerate(self.filenames):
            assert filename.endswith('.xvg')
            assert os.path.exists(path + filename), f"File {path + filename} does not exist"
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

    def get_polars(self, index):
        assert self.indexIsRetrievable(index)
        if isinstance(index, str):
            return polars.DataFrame(data=self.data[self.indexOfFilename(index)].T,
                                    schema=self.metadata[self.indexOfFilename(index)]['columns'])
        else:
            return polars.DataFrame(data=self.data[index],
                                    schema=self.metadata[index]['columns'])

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
