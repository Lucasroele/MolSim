# test_script.py
import subprocess
import filecmp
from pathlib import Path


def test_script_generates_expected_file(data_dir):
    # Paths
    input_file = data_dir / "top_cyclr/nocyc.top"
    expected_file = data_dir / "top_cyclr/expec.top"
    output_file = data_dir/ "top_cyclr/replace_me.top" # temporary file for script output

    # Run your script as a subprocess
    # Replace 'python' and 'script.py' with the correct command if needed
    subprocess.run(
        ["python3", "src/molsim/commands/top_cyclr.py", str(input_file), '-f', '-o', str(output_file)],  # add other args if needed
        check=True
    )

    # Compare the generated file to the reference
    assert filecmp.cmp(output_file, expected_file, shallow=False), "Generated file does not match expected file"

