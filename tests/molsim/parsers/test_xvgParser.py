import matplotlib.pyplot as plt
import numpy as np

from molsim.parsers.xvgParser import XvgParser
from molsim.utils.files import getFileNames


def test_XvgParser(data_dir):
    data_1 = XvgParser("xvgs/dist_1.xvg", path=data_dir)
    filepaths = []
    for f in getFileNames(".xvg", path=data_dir / "xvgs"):
        filepaths.append(data_dir / "xvgs" / f)
    assert data_1.filenames == ["xvgs/dist_1.xvg"]

    data_2 = XvgParser(filepaths)
    assert data_2.filenames[0] == "dist_1.xvg"
    assert len(data_2) == len(filepaths)
    
    first_times = []
    for dataset in data_2:
        first_times.append(dataset[0][0])
    assert first_times == [0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
    



def test_plot_scaffold_accepts_numpy_arrays_and_xvg_files(data_dir, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    values = np.array([[0, 1, 2], [0.1, 0.3, 0.2]])

    _, array_ax = XvgParser.plot_scaffold(values, xlabel="x", ylabel="y", show=False)
    assert array_ax.lines[0].get_ydata().tolist() == [0.1, 0.3, 0.2]

    filename = data_dir / "xvgs/dist_1.xvg"
    _, file_ax = XvgParser.plot_scaffold(filename, show=False)
    assert len(file_ax.lines) == 1
    assert file_ax.get_xlabel()


def test_msd_plot_uses_shared_scaffold(tmp_path, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    filename = tmp_path / "msd.xvg"
    filename.write_text(
        '@ title "Mean square displacement"\n'
        '@ xaxis label "Time (ps)"\n'
        '@ yaxis label "MSD (nm^2)"\n'
        '@ s0 legend "MSD D 1.0"\n'
        "1000 0.1\n"
        "2000 0.2\n"
    )
    parser = XvgParser(filename)

    _, ax = parser.msd_plot("test")

    assert ax.lines[0].get_xdata().tolist() == [1.0, 2.0]
    assert ax.get_xlabel() == r"$\tau$ (ns)"
    assert ax.get_ylabel() == r"MSD (nm$^2$)"


def test_rmsd_plot_uses_shared_scaffold(tmp_path, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    filenames = [tmp_path / "rmsd_1.xvg", tmp_path / "rmsd_2.xvg"]
    for filename, value in zip(filenames, (0.1, 0.3)):
        filename.write_text(
            '@ title "RMSD"\n'
            '@ xaxis label "Time (ps)"\n'
            '@ yaxis label "RMSD (nm)"\n'
            + "".join(f"{time} {value}\n" for time in (0, 1000, 2000, 3000))
        )
    parser = XvgParser(filenames)

    fig, ax = parser.rmsd_plot(
        "ensemble",
        decimals_on_xaxis=1,
        decimals_on_yaxis=3,
    )

    assert len(ax.lines) == 2
    assert ax.lines[0].get_xdata().tolist() == [0.0, 1.0, 2.0, 3.0]
    assert ax.get_xlabel() == "Time (ns)"
    assert ax.get_ylabel() == "RMSD (nm)"
    assert ax.xaxis.get_major_formatter().fmt == "%.1f"
    assert ax.yaxis.get_major_formatter().fmt == "%.3f"
    assert ax.get_xlim() == (-10.0, 1010.0)
    assert ax.get_ylim() == (0.0, 0.2)
    assert [text.get_text() for text in ax.texts] == ["0.20 ± 0.10 nm"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["ensemble"]
    plt.close(fig)


def test_rmsf_plot_shows_ensemble_mean_and_std(tmp_path, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    filenames = [tmp_path / "rmsf_1.xvg", tmp_path / "rmsf_2.xvg"]
    for filename, values in zip(filenames, ([0.1, 0.3], [0.3, 0.5])):
        filename.write_text(
            '@ title "RMS fluctuation"\n'
            '@ xaxis label "Residue"\n'
            '@ yaxis label "(nm)"\n'
            f"1 {values[0]}\n"
            f"2 {values[1]}\n"
        )
    parser = XvgParser(filenames)

    _, ax = parser.rmsf_plot("ensemble", residue_names=["A", "B"])

    assert ax.lines[0].get_ydata().tolist() == [0.2, 0.4]
    assert len(ax.collections) == 1
    assert ax.get_ylabel() == "RMSF (nm)"


def test_rmsf_plot_groups_equivalent_atoms_and_removes_hydrogens(tmp_path, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    filenames = [tmp_path / "rmsf_1.xvg", tmp_path / "rmsf_2.xvg"]
    profiles = ([0.1, 0.3, 0.2, 0.4], [0.3, 0.5, 0.4, 0.6])
    for filename, values in zip(filenames, profiles):
        filename.write_text(
            '@ title "RMS fluctuation"\n'
            '@ xaxis label "Residue"\n'
            '@ yaxis label "(nm)"\n'
            + "".join(f"{index} {value}\n" for index, value in enumerate(values))
        )
    equiv_file = tmp_path / "equivalent.json"
    equiv_file.write_text("[[0, 1], [2, 3]]")
    parser = XvgParser(filenames)

    _, ax = parser.rmsf_plot(
        "ensemble",
        equiv=equiv_file,
        residue_names=["C1", "H1"],
    )

    assert [line.get_ydata().tolist() for line in ax.lines] == [[0.2], [0.4]]
    assert len(ax.collections) == 2
    assert [label.get_text() for label in ax.get_xticklabels()] == ["C1"]
    
