import numpy as np
import pandas as pd
import polars
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
from os.path import commonprefix

class XvgParser:
    """
    Parse multiple .xvg files

    Instantiation:
        xvg = XvgParser(filenames, path=None)
    Arguments:
        filenames: str | Path | [str] | [Path]
        path: str | Path



    self.attributes
    filenames = [filenames]
    filepaths = [full paths to files]
    data      = [np.ndarray]         # Cols in rows
    metadata  = [{key, val}]
    path      = path to files
    size      = number of entries
    """

    def __init__(self, filenames, path=None, pp=True):
        self.data = []
        self.metadata = []
        # Handle path
        if isinstance(path, str):
            self.path = Path(path)
        elif isinstance(path, Path):
            self.path = path
        else:
            self.path = Path('.')
        if not self.path.is_dir():
            raise FileNotFoundError(f"The directory `{path}` does not exist.")
        # Prepare filenames
        if isinstance(filenames, str):
            self.filenames = [filenames]
            self.filepaths = [self.path / Path(filenames)]
        elif isinstance(filenames, Path):
            self.filenames = [filenames.name]
            self.filepaths = [self.path / filenames]
        elif isinstance(filenames, list) and len(filenames) > 0 and all(isinstance(f, str) for f in filenames):
            self.filenames = filenames
            self.filepaths = [self.path / Path(f) for f in filenames]
        elif isinstance(filenames, list) and len(filenames) > 0 and all(isinstance(f, Path) for f in filenames):
            self.filenames = [f.name for f in filenames]
            self.filepaths = [self.path / f for f in filenames]
        else:
            raise ValueError("First argument must be a `str`, `Path`, `[str]` or `[Path]`.")
        for f in self.filepaths:
            if f.suffix != '.xvg':
                raise TypeError("All files must have the .xvg extension.")
            if not f.is_file():
                raise FileNotFoundError(f"File `{f}` does not exist.")
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
                            if len(self.metadata[i]['columns']) + 1 == len(tempdata):
                                if self.metadata[i].get('xaxis') != None:
                                    self.metadata[i]['columns'].insert(0, self.metadata[i]['xaxis'])
                                else:
                                    self.metadata[i]['columns'].insert(0, '#')
                            elif self.metadata[i].get('xaxis') != None and \
                            self.metadata[i].get('yaxis') != None:
                                self.metadata[i]['columns'] = [
                                    self.metadata[i]['xaxis'],
                                    self.metadata[i]['yaxis']]
                                if len(self.metadata[i]['columns']) < len(tempdata):
                                    self.metadata[i]['columns'].extend(
                                        ['Extra ' + str(leg) for leg in range(1, 1 + len(tempdata) - len(self.metadata[i]['columns']))])
                            elif len(self.metadata[i]['comments']) > 0 and \
                            prev_line == self.metadata[i]['comments'][-1] and \
                            prev_line.split() == len(tempdata):
                                self.metadata[i]['columns'] = [
                                    cel for cel in prev_line.split()]
                            elif len(self.metadata[i]['columns']) == len(tempdata) - 1:
                                self.metadata[i]['columns'].insert(0, '#')
                            else:
                                self.metadata[i]['columns'] = [
                                    str(leg) for leg, _ in enumerate(tempdata)]
                        if not self.metadata[i].get('xaxis'):
                            self.metadata[i]['xaxis'] = self.metadata[i]['columns'][0]                   # Read data
                    for k, cel in enumerate(infoline):
                        tempdata[k].append(float(cel))
                    prev_line = line
            self.data.append(np.array(tempdata))
            assert len(self.metadata[i]['columns']) == len(tempdata), f"{len(self.metadata[i]['columns'])} != {len(tempdata)}"
            if pp:
                self._ppInit(i)

        self.size = len(self.data)

    def _ppInit(self, index):
        if self.metadata[index]['yaxis'] == "(nm)" and \
           self.metadata[index]['title'] == "RMS fluctuation":
            self.metadata[index]['yaxis'] = "RMSF (nm)"
            self.metadata[index]['columns'][1] = "RMSF (nm)"


    def info(self, index=None):
        def printInfo(index):
            print(f"File {index}: {self.filenames[index]}")
            print(f"  Title: {self.metadata[index].get('title', 'N/A')}")
            print(f"  X-axis: {self.metadata[index].get('xaxis', 'N/A')}")
            print(f"  Y-axis: {self.metadata[index].get('yaxis', 'N/A')}")
            print(f"  Columns: {self.metadata[index]['columns']}")
            print(f"  Number of data points: {self.data[index].shape[1]}")
            print("")

        if index is None:
            for i, filename in enumerate(self.filenames):
                printInfo(i)
        else:
            index=self._resolveIndex(index)
            printInfo(index)

    def get_numpy(self, index):
        index = self._resolveIndex(index)
        return self.data[index]

    def get_pandas(self, index):  # returns as columns of datapoints
        index = self._resolveIndex(index)
        return pd.DataFrame(data=self.data[index],
                            columns=self.metadata[index]['columns'])

    def get_polars(self, index):
        index = self._resolveIndex(index)
        return polars.DataFrame(data=self.data[index],
                                schema=self.metadata[index]['columns'])

    def plot(self, all_in_one=False, subplot_titles=None, showlegend=True, traces=None):
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
            if traces is None:
                p = len(commonprefix(self.filenames))
                s = len(commonprefix([x[::-1] for x in self.filenames]))
                traces = [x[p:max(p, len(x) - s)] for x in self.filenames]
        else:
            rows = self.size
            if subplot_titles is not None:
                if len(subplot_titles) != self.size:
                    raise ValueError(f"The argument `subplot_titles` should be of lenght: {self.size}")
            else:
                subplot_titles = [Path(f).stem for f in self.filenames]
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
                                             name=traces[i]),
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
        fig.update_layout(showlegend=showlegend)
        if not all_in_one:
            fig.update_layout(height=300 * len(self.filenames))
        fig.show()

    def rmsd_plot(self, name, cname="Blues", fname="", x_range_ns=1000, ylim=(0,0.2)):
        # def plot_individual_rmsd_xvg_files(File_path1, File_path2, File_name, File_list, File_numbers, Name, CName, Colour, Title, ylim):
        """Plot replica RMSD traces from GROMACS ``.xvg`` files and annotate summary statistics.

        For every system in ``File_list``, the RMSD time series of each replica listed in
        ``File_numbers`` is loaded and drawn on a single axes, with each replica shaded
        from a sequential colormap. Per-replica mean RMSDs and their spread are printed,
        and a summary of the mean and standard deviation is shown in the plot title. The
        figure can optionally be saved as a 300 dpi PNG.

        The path of each input file is built as::

            {File_path1}{name}{File_path2}/{File_name}{replica}.xvg

        where ``name`` is an entry of ``File_list`` and ``replica`` is a replica number
        taken from the matching entry of ``File_numbers``.

        :param File_path1: Prefix of the directory path, placed before the system name.
            Also the base directory for saved figures (``{File_path1}Figures_paper/``).
        :type File_path1: str
        :param File_path2: Suffix of the directory path, placed after the system name
            and before ``/``.
        :type File_path2: str
        :param File_name: Filename stem of the ``.xvg`` files. The replica number and
            ``.xvg`` are appended, e.g. ``"rmsd_"`` gives ``rmsd_1.xvg``.
        :type File_name: str
        :param File_list: System/directory names to process. Each is substituted for
            ``name`` in the path template above.
        :type File_list: list[str]
        :param File_numbers: Replica numbers for each system, parallel to ``File_list``.
            Each element is a comma-separated string of integers, e.g. ``"1,2,3"``.
        :type File_numbers: list[str]
        :param Name: Label shown in the legend box.
        :type Name: str
        :param CName: Name of a Matplotlib colormap (e.g. ``"Blues"``) used to colour
            individual replica traces. Sampled over the range 0.2 to 1 to avoid
            near-white colours.
        :type CName: str
        :param Colour: A single Matplotlib colour used for the title text and the legend
            line. It is independent of ``CName``, so choose a colour that matches the
            colormap.
        :type Colour: str or tuple
        :param Title: Output filename stem. If truthy, the figure is saved to
            ``{File_path1}Figures_paper/{Title}.png`` at 300 dpi. If empty or ``None``,
            nothing is saved.
        :type Title: str or None
        :param ylim: Y-axis limits as ``(ymin, ymax)`` in nm.
        :type ylim: tuple[float, float]

        :returns: ``None``. The figure is displayed with ``plt.show()`` and optionally
            written to disk.
        :rtype: None

        .. note::
        **Statistics.** Two different quantities are reported:

        * ``av_mean`` / ``av_sd``: the mean and standard deviation of all RMSD values
            pooled over every replica of a system (all frames of all replicas).
        * ``average_mean`` / ``average_std``: the mean and standard deviation of the
            *per-replica* mean RMSDs, i.e. the between-replica spread.

        The plot title shows ``av_mean ± average_std``, which mixes the pooled mean
        with the between-replica standard deviation.

        .. warning::
        Only the **last** entry of ``File_list`` contributes to the title text and the
        final printed summary, because the statistics variables are overwritten on each
        iteration of the outer loop. All systems are drawn on the same axes, but only
        the last is annotated. In practice this function is intended to be called with
        a single-system ``File_list``.

        .. warning::
        All replicas of a system must have the same number of frames, since they are
        stacked into a single array for the per-frame mean and standard deviation.
        Those per-frame quantities are currently computed but not plotted (the
        corresponding ``ax.plot`` and ``fill_between`` calls are commented out).

        **Example**::

            plot_individual_rmsd_xvg_files(
                File_path1="/data/sims/", File_path2="/analysis", File_name="rmsd_",
                File_list=["CNR4A_iBu"], File_numbers=["1,2,3,4,5"],
                Name="CNR4A + iBu", CName="Blues", Colour="tab:blue",
                Title="CNR4A_iBu_rmsd", ylim=(0.0, 4.0),
            )
        """
        def block_average(x, y, n=5):
            m = len(x) // n * n                      # drop leftover points that don't fill a block
            x_new = x[:m].reshape(-1, n).mean(axis=1)
            y_new = y[:m].reshape(-1, n).mean(axis=1)
            return x_new, y_new

        match cname:
            case "Blues":
                colour = "xkcd:vivid blue"
            case "Oranges":
                colour = "xkcd:bright orange"
            case "Greens":
                colour = "xkcd:green"
            case "Purples":
                colour = "xkcd:vivid purple"
            case "Greys":
                colour = "XKCD:black"
            case _:
                colour = "black"

        # def plot_individual_rmsd_xvg_files(File_path1, File_path2, File_name, File_list, File_numbers, Name, CName, Colour, Title, ylim):
        fig, ax = plt.subplots(figsize=(9,5))

        #for name, simulations in zip(File_list, File_numbers):
        all_variables = []
        mean_list = []
        #sim_numbers = [int(sim_num) for sim_num in simulations.split(',')]
        cmap = mpl.colormaps[cname]
        color = iter(cmap(np.linspace(0.2, 1, self.size)))

        for filename in self.filenames:
            data = self.get_numpy(filename)
            n = int(len(data[0])/500)
            x_or = (data[0]/1000)[:,None]
            y_or = data[1][:,None]
            x, y = block_average(x_or, y_or, n=n)
            all_variables.append(y_or)
            c = next(color)
            ax.plot(x, y, color = c)

            mean = np.mean(y)
            mean_list.append(mean)

        average_std = np.std(mean_list)
        average_mean = np.mean(mean_list)

        # Calculate mean and standard deviation for each time point across simulations
        concatenate = np.stack(all_variables, axis=1)
        mean_data = concatenate.mean(axis=1)
        sd_data = concatenate.std(axis=1)
        x = np.arange(len(mean_data))
        y = mean_data
        error = sd_data/2
        #ax.plot(x,y,label=f'{Name}', color=Colour[j])
        #ax.fill_between(x, y-error, y+error, alpha=0.3, color=Colour)


        ax.set_xlabel('Time (ns)',fontsize=20)
        ax.set_ylabel('RMSD (nm)',fontsize=20)
        #plt.suptitle(name, fontsize= 15)# Add subtitle above the title
        #plt.title(f'{average_mean:.2f} ± {average_std:.2f} nm', fontsize=20, x=0.16, y=0.75, color=colour) #x=0.85, y=0.01

        ax.text(0.98, 0.03, f'{average_mean:.2f} ± {average_std:.2f} nm',
        transform=ax.transAxes, fontsize=20, color=colour,
        ha='right', va='bottom')

        ax.minorticks_on()
        ax.set_ylim(ylim)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))

        ax.tick_params(which='both', width=2, labelsize=19)
        ax.tick_params(which='major', direction='inout', length=8)
        ax.tick_params(which='minor', length=4, direction = "in")
        plt.rcParams['ytick.left'] = True
        decimals_on_yaxis=2
        decimals_on_xaxis=0
        ax.yaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals_on_yaxis}f"))
        ax.xaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals_on_xaxis}f"))
        ax.set_xlim(-10, x_range_ns + 10)
        #plt.title('Average RMSD', fontsize=18)
        #fig.legend(fontsize=11, bbox_to_anchor=(1.01,1),loc='upper left')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 1, box.height])
        #fig.legend(loc = 'upper right',bbox_to_anchor=(0.99, 0.85),fontsize=20)
        custom_line = [Line2D([0], [0], color=colour, lw=2)]
        ax.legend(custom_line, [f"{name}"], bbox_to_anchor=(0.02, 0.85), loc='lower left',
                        ncols=1, borderaxespad=0., fontsize=20)
        plt.tight_layout()
        if fname:
            plt.savefig(f"{self.path}/{fname}.png", dpi=300)

        plt.show()

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

    def _validateIndex(self, index):
        if not isinstance(index, (int, str)):
            raise TypeError(f"index must be int or str, got {type(index).__name__}")
        if isinstance(index, int):
            if not 0 <= index < self.size:
                raise IndexError(f"index {index} out of range for size {self.size}")
        elif index not in self.filenames:
            raise KeyError(f"{index} is not retrievable")

    def _resolveIndex(self, index):
        self._validateIndex(index)
        return self.filenames.index(index) if isinstance(index, str) else index
