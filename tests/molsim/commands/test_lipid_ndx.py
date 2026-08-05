import subprocess
import filecmp
import shutil
import tempfile
from pathlib import Path


def test_script_generates_expected_file(data_dir):
    # Paths
    input_file = data_dir / "lipid_ndx/pull.tpr"
    expected_file_1 = data_dir/ "lipid_ndx/outsjab.ndx"
    expected_file_2 = data_dir/ "lipid_ndx/ll_added.ndx"

    tempfile_name = "temp{}.ndx"
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir) / tempfile_name.format(1)
        subprocess.run(
            ["python3",
             "src/molsim/commands/lipid_ndx.py",
             str(input_file),
             '-o', str(tmp)],  # add other args if needed
            check=True
        )
        assert filecmp.cmp(tmp, expected_file_1, shallow=False), "Generated file does not match expected file"
        tmp = Path(tmpdir) / tempfile_name.format(2)
        subprocess.run(
            ["python3",
             "src/molsim/commands/lipid_ndx.py",
             str(input_file),
             '-o', str(tmp)],  # add other args if needed
            check=True
        )
        subprocess.run(
            ["python3",
             "src/molsim/commands/lipid_ndx.py",
             str(input_file),
             '-a',
             '-ll',
             '-o', str(tmp)],  # add other args if needed
            check=True
        )
        assert filecmp.cmp(tmp, expected_file_2, shallow=False), "Generated file does not match expected file"
