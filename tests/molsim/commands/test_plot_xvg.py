import subprocess
import filecmp


def test_plot_xvg(data_dir):
    # Paths
    data_1 = [data_dir / "xvgs/dist_1.xvg"]
    data_2 = [data_dir / "xvgs/dist_2.xvg", data_dir / "xvgs/dist_3.xvg"]
    data = [data_1, data_2]
    # Run your script as a subprocess
    # Replace 'python' and 'script.py' with the correct command if needed
    for datum in data:
        subprocess.run(
                ["python3", "src/molsim/commands/plot_xvg.py"] + datum,
            check=True
        )
    subprocess.run(
            ["python3", "src/molsim/commands/plot_xvg.py", '--dir', str(data_dir / "xvgs")],
        check=True
    )

