import subprocess
import filecmp
#from pathlib import Path


def test_script_generates_expected_file(data_dir):
    # Paths
    expected_file = data_dir / "gen_cpep/AAAAA_wanted.pdb"
    output_file = data_dir/ "gen_cpep/AAAAA_cycl.pdb" # temporary file for script output

    # Run your script as a subprocess
    # Replace 'python' and 'script.py' with the correct command if needed
    subprocess.run(
        ["python3", "src/molsim/dev/gen_cpep.py", "AAAAA", '-o', str(output_file)],  # add other args if needed
        check=True
    )

    # Compare the generated file to the reference
    assert filecmp.cmp(output_file, expected_file, shallow=False), "Generated file does not match expected file"

