import numpy as np
import pandas as pd
import polars
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

class XvgParser:
    """
    Instantiation:
        xvg = XvgParser(filenames, path=None)
    Arguments:
        filenames: str | [str] | [Path]



    self.attributes
    filenames = [filenames]
    filepaths = [full paths to files]
    data      = [np.ndarray]         # Cols in rows
    metadata  = [{key, val}]
    path      = path to files
    size      = number of entries
    """

    def __init__(self, filenames, path=None):
        self.data = []
        self.metadata = []
        # Handle path
        if isinstance(path, str):
            self.path = Path(path)
        elif isinstance(path, Path):
            self.path = path
        else:
            self.path = Path('.')
        assert self.path.is_dir(), f"The directory `{path}` does not exist."
        # Prepare filenames
        if isinstance(filenames, str):
            self.filenames = [filenames]
            self.filepaths = [self.path / Path(filenames)]
        elif isinstance(filenames, list) and len(filenames) > 0 and all(isinstance(f, str) for f in filenames):
            self.filenames = filenames
            self.filepaths = [self.path / Path(f) for f in filenames]
        elif isinstance(filenames, list) and len(filenames) > 0 and all(isinstance(f, Path) for f in filenames):
            self.filenames = [f.name for f in filenames]
            self.filepaths = [self.path / f for f in filenames]
        else:
            raise ValueError("First argument must be a string or a list of strings or list of Paths.")
        for f in self.filepaths:
            assert f.suffix == '.xvg', "All files must have the .xvg extension."
            assert f.is_file(), f"File `{f}` does not exist."
        # Parse files
        for i, filename in enumerate(self.filepaths):
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
                    if line.startswith('@'):  # parse in metadata
                        infoline = line.split()
                        if len(infoline) <= 2:
                            prev_line = line
                            continue
                        # infolinen will have at least 3 elements
                        if infoline[1] == 'title':
                            self.metadata[i]['title'] = self._findString(
                                infoline[2:])
                        elif infoline[1] == 'xaxis':
                            if infoline[2] == 'label' and len(infoline) > 3:
                                self.metadata[i]['xaxis'] = self._findString(
                                    infoline[3:])
                            else:
                                self.metadata[i]['xaxis'] = self._findString(
                                    infoline[2:])
                        elif infoline[1] == 'yaxis':
                            if infoline[2] == 'label' and len(infoline) > 3:
                                self.metadata[i]['yaxis'] = self._findString(
                                    infoline[3:])
                            else:
                                self.metadata[i]['yaxis'] = self._findString(
                                    infoline[2:])
                        elif len(infoline) >= 4 and infoline[2] == "legend" and infoline[1][0] == "s":
                            self.metadata[i]['columns'].append(
                                self._findString(infoline[3:]))
                        elif len(infoline) >= 5 and infoline[1] == "legend" and infoline[2] == "string":
                            self.metadata[i]['columns'].append(
                                self._findString(infoline[4:]))
                        continue
                    # Reading data lines now
                    infoline = line.split()
                    if not reading_data:  # entered only once
                        reading_data = True
                        for cel in infoline:
                            tempdata.append([])
                        if len(self.metadata[i]['columns']) != len(tempdata):
                            if len(self.metadata[i]['comments']) > 0 and \
                            prev_line == self.metadata[i]['comments'][-1] and \
                            prev_line.split() == len(tempdata):
                                self.metadata[i]['columns'] = [
                                    cel for cel in prev_line.split()]
                            elif len(self.metadata[i]['columns']) == len(tempdata) - 1:
                                self.metadata[i]['columns'].insert(0, '#')
                            else:
                                self.metadata[i]['columns'] = [
                                    str(leg) for leg, _ in enumerate(tempdata)]
                    for k, cel in enumerate(infoline):
                        tempdata[k].append(float(cel))
                    prev_line = line
                self.data.append(np.array(tempdata))
        self.size = len(self.data)

    def info(self):
        for i, filename in enumerate(self.filenames):
            print(f"File {i}: {filename}")
            print(f"  Title: {self.metadata[i].get('title', 'N/A')}")
            print(f"  X-axis: {self.metadata[i].get('xaxis', 'N/A')}")
            print(f"  Y-axis: {self.metadata[i].get('yaxis', 'N/A')}")
            print(f"  Columns: {self.metadata[i]['columns']}")
            print(f"  Number of data points: {self.data[i].shape[1]}")
            print("")



    def get_numpy(self, index):  # returns as rows of datapoints
        assert self._indexIsRetrievable(index)
        if isinstance(index, str):
            return self.data[self._indexOfFilename(index)]
        else:  # Must be int
            return self.data[index]

    def get_pandas(self, index):  # returns as columns of datapoints
        assert self._indexIsRetrievable(index)
        if isinstance(index, str):
            print(self.metadata[self._indexOfFilename(index)]['columns'])
            return pd.DataFrame(data=self.data[self._indexOfFilename(index)].T,
                                columns=self.metadata[self._indexOfFilename(index)]['columns'])
        else:  # Must be int
            return pd.DataFrame(data=self.data[index],
                                columns=self.metadata[index]['columns'])

    def get_polars(self, index):
        assert self._indexIsRetrievable(index)
        if isinstance(index, str):
            return polars.DataFrame(data=self.data[self._indexOfFilename(index)].T,
                                    schema=self.metadata[self._indexOfFilename(index)]['columns'])
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

    def __len__(self):
        return self.size

    def __getitem__(self, index):
        return self.get_numpy(index)

    def __iter__(self):
        for i in range(self.size):
            yield self.get_numpy(i)

    def _findString(self, list_of_strings):
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

    def _indexOfFilename(self, filename):
        for i, f in enumerate(self.filenames):
            if filename == f:
                return i

    def _indexIsRetrievable(self, index):
        if isinstance(index, int):
            assert index in range(self.size)
        elif isinstance(index, str):
            assert index in self.filenames
        else:
            raise ValueError(f"{index} is not retrievable")
        return True
