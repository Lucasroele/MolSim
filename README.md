# MolSim
scripts for use preparing molecular simulations

## Installation
To install MolSim, clone the repository and navigate to the project directory:

Then, install using:

```bash
pip install -e .
```

## Usage
```bash
molsim xvg_min src/molsim/operands/dist.xvg
molsim md_check src/molsim/operands/align.gro
molsim make_posres src/molsim/operands/align.gro
molsim plot_xvg -d src/molsim/operands/xvgs
molsim frame_times src/molsim/operands/dist.xvg -b 3.5 -dr 0.5 -nf 5
molsim top_cyclr tests/fixtures/top_cyclr/nocyc.top
molsim gen_cpep AAAAA
molsim count_mols tests/fixtures/count_mols/mols.gro tests/fixtures/count_mols/before3.top -o topol.top -am "Protein_chain_A 1"
molsim split_mem src/molsim/operands/pull.tpr
```

