import numpy as np
from pathlib import Path
import re


class NdxParser:
    """
    Instantiation:
        xvg = XvgParser(filenames, path=None)
    Arguments:
        filenames: str | Path | [str] | [Path]
        path: str | Path



    self.attributes
    filenames = [filenames]
    filepaths = [full paths to files]
    path      = path of directory with files
    data      = [{group_name: np.ndarray}]
    """
    def __init__(self, filenames, path=None):
        self.data = []
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
            raise ValueError("First argument must be a `str`, `Path`, `[str]` or `[Path]`.")
        for f in self.filepaths:
            assert f.suffix == '.ndx', "All files must have the .ndx extension."
            assert f.is_file(), f"File `{f}` does not exist."
        # Start reading file
        self.data = []
        for f in self.filepaths:
            self.data.append({})
            with open(f, 'r') as file:
                sel_name = ''
                for line in file:
                    # matches [ <name> ] in the line
                    line_re = re.compile(r'^\s*\[\s*([^\]\r\n]+?)\s*\]\s*$')
                    m = line_re.match(line)
                    if line.strip() == '':
                        continue
                    if m:
                        sel_name = line.rstrip('\n').strip().lstrip('[').rstrip(']').strip()
                        self.data[-1][sel_name] = []
                        continue
                    if sel_name:
                        self.data[-1][sel_name].append(map(int, line.rstrip('\n').split()))

    def lines(self, index: int | str):
        assert self._indexIsRetrievable(index)
        if isinstance(index, str):
            index = self._indexOfFilename(index)
        for key, val in self.data[index].items():
            yield f'[ {key} ]\n'
            for row in val:
                yield ' '.join(str(i) for i in row) + '\n'
            yield '\n'

    def get_numpy(self, index: int | str):
        assert self._indexIsRetrievable(index)
        if isinstance(index, str):
            return self.data[self._indexOfFilename(index)]
        else:  # Must be int
            return self.data[index]

    def _indexOfFilename(self, filename):
        for i, f in enumerate(self.filenames):
            if filename == f:
                return i

    def _indexIsRetrievable(self, index):
        if isinstance(index, int):
            assert index in range(len(self))
        elif isinstance(index, str):
            assert index in self.filenames
        else:
            raise ValueError(f"{index} is not retrievable")
        return True

    def __len__(self):
        return len(self.filenames)

    def __getitem__(self, index):
        return self.get_numpy(index)

    def __iter__(self):
        for i in range(self.size):
            yield self.get_numpy(i)

class XvgParser:
    """
    Instantiation:
        xvg = XvgParser(filenames, path=None)
    Arguments:
        filenames: str | Path | [str] | [Path]



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
