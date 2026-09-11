import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np


class ArgumentValidationError(Exception):
    """Raised when CLI arguments fail validation."""
    pass


def register(subparsers):
    parser = subparsers.add_parser(
        'scale_coords',
        help='Scale the first coordinate file to match the box dimensions of the second coordinate file.'
    )
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(
        prog='scale_coords.py',
        description='Scale a coordinate file to the box dimensions of another coordinate file.'
    )
    parser = addArguments(parser)
    return parser.parse_args()


def addArguments(parser):
    parser.add_argument(
        'coordfile1',
        help='coordinate file to scale; trajectories use the first frame'
    )
    parser.add_argument(
        'coordfile2',
        help='coordinate file whose box dimensions define the target size'
    )
    parser.add_argument(
        '-o',
        '--output',
        default='scaled.pdb',
        help='output file name [default: scaled.pdb]'
    )
    parser.add_argument(
        '-f',
        '--force',
        action='store_true',
        help='overwrite the output file if it already exists'
    )
    return parser


def validateArguments(args):
    for attr in ('coordfile1', 'coordfile2'):
        path = Path(getattr(args, attr))
        if not path.is_file():
            raise ArgumentValidationError(f"The file `{path}` does not exist.")

    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        raise ArgumentValidationError(f"The file `{args.output}` already exists and `-f` was not used.")
    return args


def _first_frame_universe(path):
    universe = mda.Universe(str(path))
    if len(universe.trajectory) == 0:
        raise ValueError(f'No frames were found in `{path}`')
    universe.trajectory[0]
    return universe


def _validate_box_geometry(ref_u, target_u):
    ref_angles = np.asarray(ref_u.dimensions[3:6], dtype=float)
    target_angles = np.asarray(target_u.dimensions[3:6], dtype=float)

    if not np.allclose(ref_angles, target_angles):
        raise ArgumentValidationError(
            "The input files do not have the same box geometry. "
        )

    if not np.allclose(ref_angles, [90.0, 90.0, 90.0]):
        raise ArgumentValidationError(
            "The box geometry is not a rectangular box (alpha=beta=gamma=90°). "
            "This command currently only supports orthorhombic boxes; "
            "other geometries may be implemented later."
        )


def _scale_by_molecule_com(universe, scale, center):
    positions = universe.atoms.positions.copy()
    for i, residue in enumerate(universe.residues):
        atom_indices = residue.atoms.indices
        if len(atom_indices) == 0:
            continue
        mol_com = residue.atoms.center_of_mass()
        relative = universe.atoms.positions[atom_indices] - mol_com
        #scaled_com = (mol_com - center) * scale + center
        scaled_com = mol_com * scale
        positions[atom_indices] = relative + scaled_com

    if len(universe.residues) == 0:
        mol_com = universe.atoms.center_of_mass()
        relative = universe.atoms.positions - mol_com
        scaled_com = (mol_com - center) * scale + center
        positions[:] = relative + scaled_com

    return positions


def main(args):
    try:
        args = validateArguments(args)

        ref_u = _first_frame_universe(args.coordfile1)
        target_u = _first_frame_universe(args.coordfile2)
        _validate_box_geometry(ref_u, target_u)

        ref_dims = np.asarray(ref_u.dimensions[:3], dtype=float)
        target_dims = np.asarray(target_u.dimensions[:3], dtype=float)

        if np.any(ref_dims <= 0):
            raise ValueError(f'Invalid box dimensions in `{args.coordfile1}`: {ref_dims}')
        if np.any(target_dims <= 0):
            raise ValueError(f'Invalid box dimensions in `{args.coordfile2}`: {target_dims}')

        scale = target_dims / ref_dims
        center = ref_u.atoms.center_of_mass()
        scaled_positions = _scale_by_molecule_com(ref_u, scale, center)

        ref_u.atoms.positions = scaled_positions
        ref_u.dimensions = np.array([*target_dims, 90.0, 90.0, 90.0], dtype=float)
        ref_u.atoms.write(args.output)

        return args.output
    except ArgumentValidationError as exc:
        raise SystemExit(str(exc)) from exc
