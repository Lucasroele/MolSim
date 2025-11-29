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
cd src/molsim/operands | python ../plotxvg.py
```

