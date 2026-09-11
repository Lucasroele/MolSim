#!/usr/bin/env python3

import MDAnalysis as mda
import numpy as np
import sys


def scale_reference(reference, structure, output):
    ref = mda.Universe(reference)
    new = mda.Universe(structure)

    ref_box = ref.dimensions[:3]
    new_box = new.dimensions[:3]

    ref_com = ref.atoms.center_of_mass()

    # Box scaling factors
    scale = new_box / ref_box

    # Scale the reference COM with the box
    new_com = ref_com * scale

    # Keep the internal reference structure unchanged
    positions = ref.atoms.positions - ref_com + new_com

    ref.atoms.positions = positions
    ref.dimensions = new.dimensions

    ref.atoms.write(output)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <reference.gro> <structure.gro> <output.gro>")
        sys.exit(1)

    scale_reference(sys.argv[1], sys.argv[2], sys.argv[3])