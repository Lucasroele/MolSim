import subprocess
import filecmp
import tempfile
from pathlib import Path


def test_script_generates_expected_file(data_dir):
    # Paths
    input_file = data_dir / "ndxs/pull.tpr"
    expected_file = data_dir/ "ndxs/index1.ndx" # temporary file for script output

    tempfile_name = "temp.ndx"
    with tempfile.TemporaryDirectory() as tmpdir:
        output_file = Path(tmpdir) / tempfile_name

        subprocess.run(
            ["python3", "src/molsim/dev/makeLipidNDX.py", str(input_file), '-o', str(output_file)],  # add other args if needed
            check=True
        )

        # Compare the generated file to the reference
        assert filecmp.cmp(output_file, expected_file, shallow=False), "Generated file does not match expected file"

    # Run your script as a subprocess
