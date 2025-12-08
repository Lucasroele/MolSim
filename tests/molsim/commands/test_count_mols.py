# test_script.py
import subprocess
import filecmp
from pathlib import Path


def test_script_generates_expected_file(data_dir: Path):
    # Paths
    output_file = data_dir / "count_mols/replace_me.top" # temporary file for script output
    coord_file = data_dir / "count_mols/mols.gro"
    # Run your script as a subprocess
    # Replace 'python' and 'script.py' with the correct command if needed
    for i in [1,2,3]:
        input_file = data_dir / f"count_mols/before{i}.top"
        input_list = ["python3", "src/molsim/commands/count_mols.py", str(coord_file),  str(input_file), '-f', '-o', str(output_file)]
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

