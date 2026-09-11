import subprocess

import MDAnalysis as mda
import numpy as np


FIXTURE_DIR = "tests/fixtures/scale_coords"


def test_scale_coords_scales_to_target_box(tmp_path):
    ref_file = f"{FIXTURE_DIR}/to_be_scaled.gro"
    target_file = f"{FIXTURE_DIR}/new_scale_coord.gro"
    output_file = tmp_path / "scaled_from_gro.pdb"

    subprocess.run(
        ["molsim", "scale_coords", ref_file, target_file, "-o", str(output_file)],
        check=True,
    )

    output_universe = mda.Universe(str(output_file))
    target_universe = mda.Universe(target_file)
    assert np.allclose(output_universe.dimensions[:3], target_universe.dimensions[:3])
    assert output_universe.atoms.n_atoms == 21617
    assert np.isfinite(output_universe.atoms.positions).all()


def test_scale_coords_uses_first_frame_of_trajectory(tmp_path):
    ref_file = f"{FIXTURE_DIR}/to_be_scaled.gro"
    target_file = f"{FIXTURE_DIR}/new_scale_traj.pdb"
    output_file = tmp_path / "scaled_from_traj.pdb"

    subprocess.run(
        ["molsim", "scale_coords", ref_file, target_file, "-o", str(output_file)],
        check=True,
    )

    output_universe = mda.Universe(str(output_file))
    target_universe = mda.Universe(target_file)
    first_frame_dims = target_universe.trajectory[0].dimensions[:3]
    assert np.allclose(output_universe.dimensions[:3], first_frame_dims)
    assert np.isfinite(output_universe.atoms.positions).all()


def test_scale_coords_rejects_nonmatching_box_geometry(tmp_path):
    ref_file = f"{FIXTURE_DIR}/to_be_scaled.gro"
    target_file = f"{FIXTURE_DIR}/new_scale_coord.gro"
    other_target = tmp_path / "wrong_geometry.pdb"

    u = mda.Universe(str(target_file))
    u.dimensions = np.array([*u.dimensions[:3], 90.0, 80.0, 90.0], dtype=float)
    u.atoms.write(str(other_target))

    result = subprocess.run(
        ["molsim", "scale_coords", ref_file, str(other_target)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "same box geometry" in result.stderr.lower() or "same box geometry" in result.stdout.lower()
