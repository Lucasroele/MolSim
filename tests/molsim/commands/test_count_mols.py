# test_script.py
import subprocess
import filecmp
from pathlib import Path
import sys


def test_script_generates_expected_file(data_dir: Path, tmp_path: Path, pytestconfig):
    # Paths
    script = pytestconfig.rootpath / "src/molsim/commands/count_mols.py"
    coord_file = data_dir / "count_mols/mols.gro"
    for i in [1,2,3]:
        output_file = tmp_path / f"output{i}.top" # temporary file for script output
        input_file = data_dir / f"count_mols/before{i}.top"
        input_list = [sys.executable, script, str(coord_file),  str(input_file), '-f', '-o', str(output_file)]
        if i == 1:
            input_list += ['-sl', '1']
        elif i == 2:
            input_list += ['-am', 'Protein_chain_A 1 SOL 10']
        elif i == 3:
            input_list += ['-am', 'Protein_chain_A 1']
        subprocess.run(
            input_list,
            check=True
        )
        # Compare the generated file to the reference
        expected_file = data_dir / f"count_mols/after{i}.top"
        assert filecmp.cmp(output_file, expected_file, shallow=False), "Generated file does not match expected file"
    subprocess.run(
        [
            sys.executable,
            script,
            str(coord_file),
             # temporary file for script output
            '-o',
            tmp_path / "output.top"
        ],
        check=True,
        cwd=data_dir / "count_mols"
        )

