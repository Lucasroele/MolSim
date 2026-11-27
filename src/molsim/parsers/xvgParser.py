import json
from matplotlib.pylab import trace
import numpy as np
import pandas as pd
import polars
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FuncFormatter
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
    Usage:
        ```
        for d in xvg:
            print(d)  # np.ndarray of data`
        ```


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
        self.colorsets = {}
        self.colorsets["Blues"] = "xkcd:vivid blue"
        self.colorsets["Oranges"] = "xkcd:bright orange"
        self.colorsets["Greens"] = "xkcd:green"
        self.colorsets["Purples"] = "xkcd:vivid purple"
        self.colorsets["Greys"] = "xkcd:black"

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
        """
        Data is stored as a numpy.ndarray so this returns a view.
        """
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
        max_j = 0
        for i, data in enumerate(self.data):  # subplots
            for j, row in enumerate(data):  # rows
                if j == 0:
                    continue
                if j > max_j:
                    max_j = j
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
                        if trace['name'] in [t['name'] for t in fig.data]:
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

        names_colors = {}
        alt_col_counter = 0
        for t in fig.data:
            if t['name'] not in names_colors:
                cur_colors = list(dict.fromkeys(color for color in names_colors.values()))
                names_colors[t['name']] = t['line']['color']
                if t['line']['color'] in cur_colors:
                    alt_col_counter += 1
                    t['line']['color'] = colors[alt_col_counter + max_j]

        # Update layout
        fig.update_layout(showlegend=showlegend)
        # Create list of list of legend entries for each subplot
        if not all_in_one or rows > 1:
            i_height = max(300, len(names_colors) * 50 + 20)
            fig.update_layout(height=i_height * len(self.filenames))
        else:
            fig.update_layout(height=600)
        fig.show()

    @classmethod
    def plot_scaffold(
        cls,
        data,
        name="",
        cname="Blues",
        title="",
        fontsize=20,
        fname="",
        path=None,
        xlabel=None,
        ylabel=None,
        labels=None,
        stds=None,
        band_scale=0.5,
        x_limits=None,
        y_limits=None,
        tick_labels=None,
        block_size=None,
        tick_spacing=1,
        decimals_on_xaxis=None,
        decimals_on_yaxis=2,
        annotation=None,
        show=True,
    ):
        """
        This classmethod is currently being used by:
            - rmsd_plot()
            - rmsf_plot()
            - msd_plot()
        
        Plot XVG files or channel-first NumPy arrays with shared styling.

        ``data`` may be an XVG filename, a sequence of XVG filenames, a
        ``(channels, points)`` array, or a sequence of such arrays. The first
        two channels are used as x and y values. When ``stds`` is supplied,
        each corresponding series also gets a shaded standard-deviation band.
        """
        #####
        # Input verification
        #####
        parser = None
        if isinstance(data, (str, Path)) or (
            isinstance(data, (list, tuple))
            and data
            and all(isinstance(item, (str, Path)) for item in data)
        ):
            parser = cls(data, path=path)
            series = parser.data
            output_path = parser.path
            if labels is None:
                labels = [Path(filename).stem for filename in parser.filenames]
        else:
            output_path = Path(path) if path is not None else Path(".")
            try:
                if (
                    isinstance(data, (list, tuple))
                    and data
                    and all(np.asarray(values).ndim == 2 for values in data)
                ):
                    series = [np.asarray(values, dtype=float) for values in data]
                else:
                    array = np.asarray(data, dtype=float)
                    if array.ndim == 2:
                        series = [array]
                    elif array.ndim == 3:
                        series = list(array)
                    else:
                        raise ValueError
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "data must be XVG filenames or numeric NumPy arrays."
                ) from error

        if not series:
            raise ValueError("At least one data series is required.")
        for index, values in enumerate(series):
            if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] == 0:
                raise ValueError(
                    f"Series {index} must have at least two channels and one data point."
                )
        #####
        # Label handling
        #####
        single_label = False
        if labels is None:
            labels = (
                [name]
                if len(series) == 1 and name
                else [f"Run {index + 1}" for index in range(len(series))]
            )
        elif isinstance(labels, str):
            single_label = labels
        elif len(labels) != len(series):
            raise ValueError(f"Got {len(labels)} labels for {len(series)} data series.")

        #####
        # Standard deviation handling
        #####
        if stds is not None:
            std_array = np.asarray(stds, dtype=float)
            if std_array.ndim == 1 and len(series) == 1:
                std_series = [std_array]
            elif std_array.ndim == 2:
                std_series = list(std_array)
            else:
                raise ValueError(
                    "stds must contain one one-dimensional array per data series."
                )
            if len(std_series) != len(series):
                raise ValueError(
                    f"Got {len(std_series)} standard-deviation arrays for "
                    f"{len(series)} data series."
                )
            for index, (values, deviations) in enumerate(zip(series, std_series)):
                if len(deviations) != values.shape[1]:
                    raise ValueError(
                        f"Series {index} has {values.shape[1]} points but "
                        f"{len(deviations)} standard deviations."
                    )
        else:
            std_series = None
        #####
        # Tick label handling
        #####
        if tick_labels is not None:
            if any(values.shape[1] != series[0].shape[1] for values in series) or (
                len(tick_labels) != series[0].shape[1]
            ):
                raise ValueError(
                    "tick_labels must match the points in each data series."
                )
            if tick_spacing < 1:
                raise ValueError("tick_spacing must be a positive integer.")

        #####
        # Block average
        #####
        def block_average(x, y, n=5):
            m = len(x) // n * n                      # drop leftover points that don't fill a block
            x_new = x[:m].reshape(-1, n).mean(axis=1)
            y_new = y[:m].reshape(-1, n).mean(axis=1)
            return x_new, y_new
        #####
        # Create and populate the figure and axes
        #####
        colors = mpl.colormaps[cname](np.linspace(0.2, 1, len(series)))
        fig, ax = plt.subplots(
            figsize=(12, 6) if tick_labels is not None else (9, 5)
        )
        seen_labels = set()
        for index, (values, label, color) in enumerate(zip(series, labels, colors)):
            if block_size is not None:
                x_values, y_values = block_average(values[0], values[1], n=block_size)
            else:
                x_values, y_values = values[0], values[1]
            display_label = label if label not in seen_labels else "_nolegend_"
            seen_labels.add(label)
            line, = ax.plot(x_values, y_values, color=color, label=display_label)
            if std_series is not None:
                deviation = std_series[index] * band_scale
                ax.fill_between(
                    x_values,
                    y_values - deviation,
                    y_values + deviation,
                    color=line.get_color(),
                    alpha=0.3,
                )

        if xlabel is None and parser is not None:
            xlabel = parser.metadata[0].get("xaxis")
        if ylabel is None and parser is not None:
            ylabel = parser.metadata[0].get("yaxis")
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if x_limits is not None:
            ax.set_xlim(x_limits)
        if y_limits is not None:
            ax.set_ylim(y_limits)
        if tick_labels is not None:
            x_values = series[0][0]
            ax.set_xticks(x_values[::tick_spacing])
            ax.set_xticklabels(
                list(tick_labels)[::tick_spacing],
                rotation=45,
                fontsize=fontsize,
            )
            if not isinstance(tick_labels[0], str):
                ax.xaxis.minorticks_on()
                ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.minorticks_on()
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(which="both", width=2, labelsize=fontsize)
        ax.tick_params(which="major", direction="inout", length=8)
        ax.tick_params(which="minor", length=4, direction="in")
        if decimals_on_xaxis is not None:
            ax.xaxis.set_major_formatter(
                FormatStrFormatter(f"%.{decimals_on_xaxis}f")
            )
        ax.yaxis.set_major_formatter(
            FormatStrFormatter(f"%.{decimals_on_yaxis}f")
        )
        #####
        # Legend creation
        #####
        if single_label:
            cmap = mpl.colormaps[cname]
            # `color` for individual RMSD traces
            color = cmap(np.linspace(0.2, 1, len(series)))
            # `colour` for annotation text and the legend line
            colour = color[len(series) // 2]
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 1, box.height])
            #fig.legend(loc = 'upper right',bbox_to_anchor=(0.99, 0.85),fontsize=20)
            custom_line = [Line2D([0], [0], color=colour, lw=2)]
            ax.legend(custom_line, [f"{single_label}"], bbox_to_anchor=(0.02, 0.85), loc='lower left',
                            ncols=1, borderaxespad=0., fontsize=20)
        else:
            ax.legend(fontsize=fontsize)
        #####
        # Add info in the bottom right corner?
        #####
        if annotation:
            ax.text(
                0.98,
                0.03,
                annotation,
                transform=ax.transAxes,
                fontsize=fontsize,
                color=colors[len(colors) // 2],
                ha="right",
                va="bottom",
            )
        if title:
            ax.set_title(title, fontsize=fontsize)
        fig.tight_layout()
        if fname:
            fig.savefig(output_path / f"{fname}.png", dpi=300)
        if show:
            plt.show()
        return fig, ax

    def rmsd_plot(
        self,
        name,
        cname="Blues",
        title="",
        fname="",
        x_range_ns=1000,
        ylim=(0, 0.2),
        fontsize=20,
        labels=None,
        decimals_on_xaxis=0,
        decimals_on_yaxis=2,
        block_size=None
    ):

        if labels == None:
            labels = [name] * self.size
        series = []
        mean_list = []
        for filename in self.filenames:
            data = self.get_numpy(filename)
            if block_size is None:
                block_size = max(1, int(len(data[0]) / 500))
            x = (data[0]/1000)[:,None]
            y = data[1][:,None]
            mean = np.mean(y)
            mean_list.append(mean)
            series.append(np.vstack((x.T, y.T)))

        average_std = np.std(mean_list)
        average_mean = np.mean(mean_list)

        return self.plot_scaffold(
            series,
            name=name,
            cname=cname,
            title=title,
            fontsize=fontsize,
            fname=fname,
            path=self.path,
            xlabel="Time (ns)",
            ylabel="RMSD (nm)",
            labels=labels,
            x_limits=(-x_range_ns / 100, x_range_ns + x_range_ns / 100),
            y_limits=ylim,
            decimals_on_xaxis=decimals_on_xaxis,
            decimals_on_yaxis=decimals_on_yaxis,
            annotation=f"{average_mean:.2f} ± {average_std:.2f} nm",
            block_size=block_size,
        )

    def rmsf_plot(
        self,
        name,
        cname="Blues",
        fname="",
        equiv=None,
        ylim=(0, 0.2),
        title="",
        fontsize=15,
        tick_labels=None,
        labels=None
    ):
        """Plot RMSF profiles, optionally averaging groups of equivalent atoms.

        ``equiv`` contains zero-based atom-index groups. Groups labeled with an atom
        name starting with ``H`` are omitted.
        """

        #####
        # Do some input verification and preprocessing
        #####
        lengths = [data.shape[1] for data in self.data]
        if any(length != lengths[0] for length in lengths[1:]):
            mismatched = next(
                i for i, length in enumerate(lengths) if length != lengths[0]
            )
            raise ValueError(
                f"All .xvg files must have the same number of data points. "
                f"File {self.filenames[mismatched]} has {lengths[mismatched]} points, "
                f"while the first file has {lengths[0]}."
            )

        x_values = self.data[0][0]
        for index, data in enumerate(self.data[1:], start=1):
            if not np.array_equal(data[0], x_values):
                raise ValueError(
                    f"File {self.filenames[index]} has x values that do not match "
                    "the first file."
                )
        #####
        # Handle equivalent atom groups if provided
        #####
        if equiv is not None:
            if isinstance(equiv, (str, Path)):
                with open(equiv) as file:
                    equiv = json.load(file)
            if not isinstance(equiv, (list, tuple)):
                raise ValueError(
                    "Argument `equiv` must be a list, tuple, or JSON filename."
                )
            # groups contains the atom indices for each equivalent atom group
            groups = []
            for group_index, group in enumerate(equiv):
                if not isinstance(group, (list, tuple, np.ndarray)):
                    raise ValueError(
                        f"Equivalent atom group {group_index} must be a sequence."
                    )
                indices = np.asarray(group)
                if indices.ndim != 1 or indices.size == 0:
                    raise ValueError(
                        f"Equivalent atom group {group_index} must not be empty."
                    )
                if not np.issubdtype(indices.dtype, np.integer):
                    raise ValueError(
                        f"Equivalent atom group {group_index} must contain integer indices."
                    )
                if np.any(indices < 0) or np.any(indices >= lengths[0]):
                    raise ValueError(
                        f"Equivalent atom group {group_index} has an index outside "
                        f"the available range 0..{lengths[0] - 1}."
                    )
                groups.append(indices.astype(int))

            if not groups:
                raise ValueError("equiv must contain at least one atom group.")

            if tick_labels is None:
                tick_labels = [f"Group {index + 1}" for index in range(len(groups))]
            elif len(tick_labels) != len(groups):
                raise ValueError(
                    f"tick_labels must contain {len(groups)} tick_labels, "
                    f"got {len(tick_labels)}."
                )

            keep_groups = [
                index
                for index, atom_name in enumerate(tick_labels)
                if not str(atom_name).startswith("H")
            ]
            if not keep_groups:
                raise ValueError("No non-hydrogen equivalent atom groups remain to plot.")

            x_plot = np.arange(len(keep_groups))
            means_list = []
            stds_list = []
            for data in self.data:
                means = np.asarray(
                    [np.mean(data[1][indices]) for indices in groups]
                )
                deviations = np.asarray(
                    [np.std(data[1][indices]) for indices in groups]
                )
                means_list.append(means[keep_groups])
                stds_list.append(deviations[keep_groups])
            series = [
                np.vstack((x_plot, means)) for means in means_list
            ]
            if labels is None:
                labels = [Path(filename).stem for filename in self.filenames]
            return self.plot_scaffold(
                series,
                name=name,
                cname=cname,
                title=title,
                fontsize=fontsize,
                fname=fname,
                path=self.path,
                xlabel="Names",
                ylabel="RMSF (nm)",
                labels=labels,
                stds=stds_list,
                band_scale=0.5,
                y_limits=ylim,
                tick_labels=[tick_labels[index] for index in keep_groups],
                show=True,
            )

        #####
        # No equivalent atom groups, just plot the RMSF profiles directly
        #####
        if tick_labels is not None and len(tick_labels) != lengths[0]:
            raise ValueError(
                f"labels must contain {lengths[0]} labels, "
                f"got {len(tick_labels)}."
            )

        y_values = np.asarray([data[1] for data in self.data])
        mean_values = np.mean(y_values, axis=0)
        std_values = np.std(y_values, axis=0)
        if labels is None:
            labels = [name]
        return self.plot_scaffold(
            np.vstack((x_values, mean_values)),
            name=name,
            cname=cname,
            title=title,
            fontsize=fontsize,
            fname=fname,
            path=self.path,
            xlabel="Residue",
            ylabel="RMSF (nm)",
            labels=labels,
            stds=std_values,
            y_limits=ylim,
            tick_labels=tick_labels,
            show=True,
        )

    def msd_plot(self, name, cname="Blues", title="", fontsize=20, fname="", x_range_ns=None, ylim=(0, 0), ylabel=''):
        """
        Plot MSD traces from GROMACS ``.xvg`` files and annotate summary statistics.
        """
        #####
        # Do some input verification and preprocessing
        #####
        if x_range_ns is None:
            x_range_ns = max([max(data[0]) for data in self.data]) / 1000
        if ylim == (0, 0):
            y_max = max([max(data[1]) for data in self.data])
            ylim = (0, y_max * 1.1)
        D_vals = [
            float(metadata["columns"][1].split()[2])
            for metadata in self.metadata
        ]
        D_mean = np.mean(D_vals)
        D_std = np.std(D_vals)
        if D_mean == 0:
            raise ValueError("The mean diffusion coefficient must be nonzero.")
        exponent = int(np.floor(np.log10(abs(D_mean))))
        scale = 10**exponent
        print_exponent = -5 + exponent
        annotation = rf'$D = ({D_mean/scale:.2f} \pm {D_std/scale:.2f})\times 10^{{{print_exponent}}}\ \mathrm{{cm}}^2$s$^{{{-1}}}$'
        series = [np.vstack((data[0] / 1000, data[1])) for data in self.data]
        return self.plot_scaffold(
            series,
            name=name,
            cname=cname,
            title=title,
            fontsize=fontsize,
            fname=fname,
            path=self.path,
            xlabel=r"$\tau$ (ns)",
            ylabel=ylabel or r"MSD (nm$^2$)",
            labels=[name] * self.size,
            x_limits=(-x_range_ns / 100, x_range_ns + x_range_ns / 100),
            y_limits=ylim,
            annotation=annotation,
        )

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


def plot_rmsf_multi(means_list, stds_list, residue_names, cname="Blues", fname="", labels=None, tick_spacing=1,
                    ax=None, band_scale=0.5, title=None, fontsize=15):
    """
    Plot several per-residue RMSF profiles (mean line + shaded SD band) on one axis.

    means_list, stds_list : lists of sequences, e.g. 5 arrays of 100 values each
    residue_names         : sequence of residue labels, shared by all runs (length 100)
    labels                : legend label per run (defaults to "Run 1", "Run 2", ...)
    tick_spacing          : label every Nth residue on the x axis
    ax                    : existing Axes to draw on; a new figure is created if None
    band_scale            : fraction of the SD used for the shaded band
    """
    if len(means_list) != len(stds_list):
        raise ValueError(f"Got {len(means_list)} mean arrays but {len(stds_list)} std arrays")

    n_res = len(residue_names)
    for i, (m, s) in enumerate(zip(means_list, stds_list)):
        if not (len(m) == len(s) == n_res):
            raise ValueError(
                f"Run {i}: means={len(m)}, stds={len(s)}, residue_names={n_res} (must all match)"
            )

    def diffSubstr(strings=[]):
        """
        Strips the right and left overlapping parts to return the middle that is differing.
        Might still return ["1c1","2c2"]
        """
        p = len(commonprefix(strings))
        s = len(commonprefix([x[::-1] for x in strings]))
        return [x[p:max(p, len(x) - s)] for x in strings]
    if labels is None:
        labels = diffSubstr([f"Run {i + 1}" for i in range(len(means_list))])
    elif len(labels) != len(means_list):
        raise ValueError(f"Got {len(labels)} labels for {len(means_list)} runs")
    else:
        labels = diffSubstr(labels)

    colorsets = ["Blues", "Oranges", "Greens", "Purples", "Greys"]
    if cname not in colorsets:
        raise ValueError(f"Color set '{cname}' not recognized. Choose from {colorsets}.")

    cmap = mpl.colormaps[cname](np.linspace(0.2, 1, len(means_list)))


    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    x_names = np.arange(n_res)
    for means, stds, lab, color in zip(means_list, stds_list, labels, cmap):
        means = np.asarray(means, dtype=float)
        error = np.asarray(stds, dtype=float) * band_scale
        line, = ax.plot(x_names, means, label=lab, color=color)
        ax.fill_between(x_names, means - error, means + error, alpha=0.3, color=line.get_color())
    ax.set_xticks(x_names[::tick_spacing])
    ax.set_xticklabels(list(residue_names)[::tick_spacing], rotation=45, fontsize=fontsize)
    ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    ax.tick_params(axis="y", labelsize=fontsize)
    ax.set_xlabel("Names", fontsize=fontsize)
    ax.set_ylabel("RMSF (nm)", fontsize=fontsize)
    ax.legend(fontsize=13)
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which='both', width=2, labelsize=19)
    ax.tick_params(which='major', direction='inout', length=8)
    ax.tick_params(which='minor', length=4, direction = "in")
    plt.rcParams['ytick.left'] = True
    ax.yaxis.set_major_formatter(FuncFormatter(
    lambda y, _: "0" if y == 0 else f"{y:.2f}"
    ))
    plt.tight_layout()
    if title:
        ax.set_title(title, fontsize=fontsize)
    if fname:
        # save inside dir called `rmsf` if it exists, otherwise save in current dir
        rmsf_dir = Path("rmsf")
        if rmsf_dir.exists():
            if rmsf_dir.is_dir():
                plt.savefig(rmsf_dir / f"{fname}.png", dpi=300)
        else:
            plt.savefig(f"{fname}.png", dpi=300)

    plt.show()