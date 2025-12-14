import subprocess
import filecmp
import shutil
import tempfile
from pathlib import Path


def test_script_generates_expected_file(data_dir):
    # Paths
    input_file = data_dir / "split_mem/pull.tpr"
    expected_file = data_dir/ "split_mem/sjab.ndx" # temporary file for script output

    tempfile_name = "temp.ndx"
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir) / tempfile_name
        for i in range(1,2):
            output_file = data_dir / f"split_mem/index{i}.ndx" 

            shutil.copy2(output_file, tmp)

            subprocess.run(
                ["python3", "src/molsim/commands/split_mem.py", str(input_file), '-f', '-a', '-o', str(output_file)],  # add other args if needed
                check=True
            )

            # Compare the generated file to the reference
            assert filecmp.cmp(output_file, expected_file, shallow=False), "Generated file does not match expected file"

            shutil.copy2(tmp, output_file)
        # Run your script as a subprocess
