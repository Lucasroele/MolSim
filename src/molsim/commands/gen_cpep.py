import Bio.PDB
import PeptideBuilder
from PeptideBuilder import Geometry
import argparse
import math
import numpy as np
from pathlib import Path

from molsim.utils.maths import findLstsqPlane
from molsim.utils.maths import rotationMatrix
from molsim.utils.maths import angleBetween
from molsim.utils.maths import distConseqPoints
from molsim.utils.maths import cartToCyl
from molsim.utils.maths import cylToCart


def register(subparsers):
    parser = subparsers.add_parser('gen_cpep',
                                   help='Creates a pdb file with a cyclic peptide oriented in space as a ring with its terminal nitrogen and oxygen overlapping.')
    parser = addArguments(parser)
    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='gen_cpep.py',
                                     description='Creates a pdb file with a cyclic peptide oriented in space as a ring with its terminal nitrogen and oxygen overlapping.',
                                     epilog='')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def addArguments(parser):
    parser.add_argument('aastring')
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        nargs='?',
                        default=None,
                        help='name of the output which defaults to `<aastring>_cycl.pdb`')
    parser.add_argument('-id',
                        '--newChainID',
                        type=str,
                        default='X',
                        help='character to put in column 22 in the new pdb file corresponding to the chainID.')
    parser.add_argument('-c',
                        '--center',
                        action='store_true',
                        help='centers everything in the box using the com of the points in space')
    parser.add_argument('-ap',
                        '--allpos',
                        action='store_true',
                        help='makes sure all atoms have only positive coordinates and shifts the system to satisfy this when necessary')
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='turn on MDAnalysis warnings')
    parser.add_argument('-f',
                        '--force',
                        action='store_true',
                        help='force overwrite of output file')
    parser.add_argument('-hn',
                        '--termHN',
                        action='store_true',
                        help='add a hydrogen to the terminal N')
    return parser


def validateArguments(args):
    assert len(args.aastring) > 4, "Cyclic peptides are physically unstable with less than 4 amino acids."
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    assert all(aa.upper() in all_aa for aa in args.aastring), "Please provide a valid string of amino acid residue letters."
    if args.output is not None:
        out_path = Path(args.output)
        assert not out_path.is_file() or args.force, "Output file already exists. Use -f to overwrite."
        if not args.output.endswith(".pdb"):
            args.output = args.output + ".pdb"
        else:
            args.output = args.output
    else:
        args.output = f"{args.aastring}_cycl.pdb"
    return args


def main(args):
    """
    Creates a pdb file with a peptide in a ring conformation.
    """
    args = validateArguments(args)


    offset = 360 / (len(args.aastring) - 0.5) * math.sqrt(3) / 2
#    if args.allAngles is not None:
#        args.phi = args.allAngles
#        args.psi_im1 = args.allAngles
#        args.omega = args.allAngles
#
    #offset = 90
    geo = Geometry.geometry(args.aastring[0])
    geo.phi = 180 + offset
    geo.psi = 180 + offset
    structure = PeptideBuilder.initialize_res(geo)
    for i, aa in enumerate(args.aastring):
        if i == 0:
            continue
        if i % 2 == 1:
            phi = 180 - offset
            psi = 180 - offset
        else:
            phi = 180 + offset
            psi = 180 + offset
        geo = Geometry.geometry(aa)
        geo.phi = phi
        geo.psi_im1 = psi
        PeptideBuilder.add_residue(structure, geo)
    #PeptideBuilder.add_terminal_OXT(structure)

    mol = next(iter(structure[0]))
    coords = np.array([a.coord for a in structure.get_atoms()])
    CCN_coords = np.array(
        [a.coord for a in structure.get_atoms() if a.id in ['N', 'C', 'CA']])

    for atom in next(iter(mol)):
        if atom.id == 'N':
            firstNcoords = np.array([atom.coord])
    cog_ring = np.average(CCN_coords, 0)

    CCN_coords -= cog_ring
    coords -= cog_ring
    firstNcoords -= cog_ring

    # Move into xy plane
    fit, order = findLstsqPlane(CCN_coords)
    vecs_in_plane = np.zeros((2, 3))  # use for normal_to_plane finding
    for i in range(2):
        vecs_in_plane[i][i] = 1
        vecs_in_plane[i][2] = fit[0] * \
            vecs_in_plane[i][0] + fit[1] * vecs_in_plane[i][1]
    backorder = np.array(range(3))
    for i, d in enumerate(order):
        backorder[d] = i
    # order indexes assignment of xyz into fit
    # backorder indexes assignment of fit into xyz
    normal_to_plane = np.cross(vecs_in_plane[0], vecs_in_plane[1])
    normal_to_plane = np.array(
        [normal_to_plane[backorder[0]], normal_to_plane[backorder[1]], normal_to_plane[backorder[2]]])
    zvec = np.array([0, 0, 1])
    axis = np.cross(normal_to_plane, zvec)
    angle = angleBetween(normal_to_plane, zvec)
    rot_mat = rotationMatrix(axis, angle)

    coords = np.matmul(coords, rot_mat)
    CCN_coords = np.matmul(CCN_coords, rot_mat)
    firstNcoords = np.matmul(firstNcoords, rot_mat)
    # in xy plane
    ###
    # Work with cyl coords from here
    ###

    cyl_coords = cartToCyl(coords)
    cyl_CCN_coords = cartToCyl(CCN_coords)
    cyl_firstNcoords = cartToCyl(firstNcoords)

    # Put N at y=0
    cyl_coords[:, 1] -= cyl_firstNcoords[:, 1]
    cyl_CCN_coords[:, 1] -= cyl_firstNcoords[:, 1]

    # Mirror everything through yz plane if second angle is smaller than first angle
    if cyl_CCN_coords[1, 1] - cyl_CCN_coords[0, 1] < 0:
        cyl_CCN_coords[:, 1] *= -1

    # Use only positive angles
    for i, vec in enumerate(cyl_CCN_coords):
        if i == 0:
            continue
        if vec[1] < 0:
            vec[1] += 2 * np.pi

    ang = cyl_CCN_coords[1:-1, 1] - cyl_CCN_coords[0:-2, 1]

    av_ang = np.average(ang)
    for a in ang:
        assert a > 0, "The script doesn'tn really make sense since the angle along the backbone is not monotonous."
        if np.abs(a - av_ang) < 1e-1:
            "The angle along the backbone has a high variance."
    av_dist = np.average(cyl_CCN_coords[:, 0])
    portions = len(ang) + 2
    scale_ang = (2 * np.pi) / (portions * av_ang)

    # For isosceles
    # b = 2 * av_dist * sin(av_ang / 2)
    new_av_dist = av_dist * np.sin(av_ang / 2) / np.sin(scale_ang * av_ang / 2)
    dist_diff = new_av_dist - av_dist
    scale_dist = new_av_dist / av_dist

    m_a = {}
    m_a["a"] = np.empty_like(cyl_coords)
    m_a["m"] = np.empty_like(cyl_coords)
    for i in range(len(cyl_coords)):
        m_a["m"][i][0] = cyl_coords[i][0] * scale_dist
        m_a["m"][i][1] = cyl_coords[i][1] * scale_ang
        m_a["m"][i][2] = cyl_coords[i][2]
        m_a["a"][i][0] = cyl_coords[i][0] + dist_diff
        m_a["a"][i][1] = cyl_coords[i][1] * scale_ang
        m_a["a"][i][2] = cyl_coords[i][2]
    dist_before = distConseqPoints(cylToCart(cyl_coords))
    m_a["m"] = cylToCart(m_a["m"])
    m_a["a"] = cylToCart(m_a["a"])
    dist_after = {}
    dist_after["m"] = distConseqPoints(m_a["m"])
    dist_after["a"] = distConseqPoints(m_a["a"])

    # Choose the better, delete the other
    if np.sum(np.abs(dist_after["m"] - dist_before)) < np.sum(np.abs(dist_after["a"] - dist_before)):
        use_m_or_a = "m"
        del m_a["a"]
    else:
        use_m_or_a = "a"
        del m_a["m"]
    for i, atom in enumerate(structure.get_atoms()):
        atom.coord = m_a[use_m_or_a][i]
    b = 2 * av_dist * np.sin(av_ang / 2)

    # Add a hydrogen to the first amino acid so it is easier to make cyclic peptide
    if args.termHN:
        for residue in structure.get_residues():
            first_res = residue
            break
        for atom in first_res.get_atoms():
            if atom.id == 'N':
                new_atom = atom.copy()
                new_atom.coord[0] += b
                new_atom.id = 'H'
                new_atom.element = 'H'
                new_atom.mass = 1.008
                new_atom.name = 'H'
                new_atom.fullname = 'H'
                break
        first_res.insert(1, new_atom)
        print("A hydrogen atom was added to the first residue in the chain labeled `HN`")
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(args.output)

    print(
        f"A 3D structure of the peptide represented by amino acid sequence `{args.aastring.upper()}` was stored in `{args.output}`.")


if __name__ == "__main__":
    args = parseArguments()
    main(args)
