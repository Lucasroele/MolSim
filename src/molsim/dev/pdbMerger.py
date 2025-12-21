import MDAnalysis as mda
from MDAnalysis.transformations import boxdimensions
import os
import sys
import argparse
import numpy as np
import math
from collections import Counter
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly
from scipy.linalg import lstsq
from numba import njit
from time import perf_counter
import time
from functools import lru_cache
import warnings
import contextlib

def getVirtBox(u):
    virtualbox = []
    for d in range(0,3):
        virtualbox.append(round(max(u.coord.positions[:,d]), 3))
        mind = round(min(u.coord.positions[:,d]), 3)
        if mind < 0:
            virtualbox[d] -= mind
    return virtualbox


def getChainIDs(u):
    ids = []
    for atom in u.atoms:
        if hasattr(atom, 'chainID') and atom.chainID not in ids:
            ids.append(atom.chainID)
    return ids


def angleBetween(v1, v2):
    """
    Returns the angle between 2 vectors
    """
    v1_u = v1 / np.linalg.norm(v1)  # normalizes v1
    v2_u = v2 / np.linalg.norm(v2)  # normalizes v2
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def rotationMatrix(axis, theta):
    """
    returns rotation matrix
    theta in radians
    Rotation follows the right hand rule with thumb along axis direction
    """
    axis = np.asarray(axis.flatten())
    axis = -axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def plotVecs(vecset, add_lines=True):
    """
    Make plots of n*3 sized np.array coordinates with lines from the origin outwards
    """
    if isinstance(vecset, dict):
        num = len(vecset)
    else:
        vecset = {"vectors": vecset}
        num = 1
    linedict = {}
    max_xyz = [0, 0, 0]
    specs = [[{'type': 'scatter3d'} for vec in vecset]]
    fig = make_subplots(rows=1, cols=num, specs=specs)
    i = 0
    trace1, trace2 = {}, {}
    for key, vecs in vecset.items():
        linedict[key] = [[], [], []]
        for d in range(3):
            curmax = max([max(vecs[:, d]), -min(vecs[:, d])])
            if curmax > max_xyz[d]:
                max_xyz[d] = curmax
        for vec in vecs:
            for d in range(3):
                linedict[key][d].append(0)
                linedict[key][d].append(vec[d])
                linedict[key][d].append(None)
        trace1[key] = go.Scatter3d(x=vecs[:, 0],
                                   y=vecs[:, 1],
                                   z=vecs[:, 2],
                                   mode='markers',
                                   name='markers')
        fig.add_trace(trace1[key], row=1, col=i + 1)
        if add_lines:
            trace2[key] = go.Scatter3d(x=linedict[key][0],
                                       y=linedict[key][1],
                                       z=linedict[key][2],
                                       mode='lines',
                                       name='lines')
            fig.add_trace(trace2[key], row=1, col=i + 1)
        i += 1

    for i in range(num):
        #fig.update_scenes(xaxis_range=[-max(max_xyz), max(max_xyz)],
        #                          yaxis_range=[-max(max_xyz), max(max_xyz)],
        #                          zaxis_range=[-max(max_xyz), max(max_xyz)],
        #                          row=1,
        #                          col=i+1)
        fig.update_scenes(xaxis_range=[-max_xyz[0], max_xyz[0]],
                          yaxis_range=[-max_xyz[1], max_xyz[1]],
                          zaxis_range=[-max_xyz[2], max_xyz[2]],
                          row=1,
                          col=i + 1)
#    fig.show()
    plotly.offline.iplot(fig, filename='simple-3d-scatter')
    return


def FindLstsqPlane(cloud=None, order=None, show_output=False, axis=0, test=False):
    """
    cloud should be a numpy.array
    order is a three element iterable where in terms of dimensions fit(order[0,1]) = order[2]
    axis only used if order not specified
    axis tells which axis should be in the plane 0=x 1=y 2=z
    show_output=True will plot the cloud with the plane and print some other usefull things
    the returned order is which columns are used at the respective index
    fit[0] * cloud[order[0]] + fit[1] * cloud[order[1]] + fit[2] ~= cloud[order[2]]
    """
    xyz = ['x','y','z']
    standard_order = [0,1,2]
    others = [1,2]
    well_defined = False
    # Only used when testing
    TARGET_X_SLOPE, TARGET_Y_SLOPE = int, int
    if cloud is None:
        show_output = True
        if not test:
            print("No cloud of points was specified, so one will be generated.")
    if test or cloud is None:
        # These constants are to create random data for the sake of this example
        N_POINTS = 10
        TARGET_X_SLOPE = 2
        TARGET_Y_SLOPE = 1
        TARGET_OFFSET  = 5
        EXTENTS = 5
        NOISE = 5
        # Create random data.
        # In your solution, you would provide your own xs, ys, and zs data.
        xs = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
        ys = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
        zs = []
        for i in range(N_POINTS):
            zs.append(xs[i]*TARGET_X_SLOPE + \
                      ys[i]*TARGET_Y_SLOPE + \
                      TARGET_OFFSET + np.random.normal(scale=NOISE))
        cloud = np.array(np.vstack((xs,ys,zs))).T
    try:
        assert cloud.shape[1] == 3
    except TypeError:
        print("The function FindLstsqPlane needs an n*3 size np.array to work with")
        return None, None
    # parse order (3rd element is the distance to the plane that should be minimized)
    if not order is None:
        well_defined = True
        try:  # can len() be called?
            if len(order) == 3:
                for d in order:
                    if not d in standard_order:
                        well_defined = False
                        break
            else:
                well_defined = False
        except TypeError:
            well_defined = False
    if not well_defined:
        order = standard_order
    # define the order based on axis
    new_order = list(order)
    if not well_defined and axis != 0:
        new_order[0] = axis
        i = 1
        for d in standard_order:
            if axis != d:
                others[i-1] = d
                i += 1
        if max(cloud[:,others[0]]) - min(cloud[:,others[0]]) < max(cloud[:,others[1]]) - min(cloud[:,others[1]]):
            new_order[1] = others[0]
            new_order[2] = others[1]
        else:
            new_order[1] = others[1]
            new_order[2] = others[0]
    order = new_order
    # What matters
    cloud = np.hstack((cloud[:,[order[0]]], cloud[:,[order[1]]], cloud[:,[order[2]]]))
    # Original columns are now at the index of order
    # [2,0,1] means original column 2 is now first
    A = np.hstack((cloud[:,[0,1]], np.ones((cloud.shape[0],1))))
    b = cloud[:,[2]]
    fit, residual, rnk, s = lstsq(A, b)

    # Show the points and plane for verification
    if show_output:
        if test:
            print("Points were generated along a plane with:")
            print(f"    x_slope = {TARGET_X_SLOPE}")
            print(f"    y_slope = {TARGET_Y_SLOPE}")
        fig = make_subplots(rows=1, cols=1, specs=[[{'type': "scatter3d"}]])
        lim_xyz = np.empty((3,2), dtype=float)
        for d in range(3):
            lim_xyz[d][0] = min(cloud[:,d])
            lim_xyz[d][1] = max(cloud[:,d])
        print("solution: %f %c + %f %c + %f = %c" % (fit[0], xyz[order[0]], fit[1], xyz[order[1]], fit[2], xyz[order[2]]))
        print("residual:", residual)
        cloud_trace = go.Scatter3d(x=cloud[:,0].flatten(),
                                   y=cloud[:,1].flatten(),
                                   z=cloud[:,2].flatten(),
                                   mode='markers',
                                   name='markers')
        fig.add_trace(cloud_trace, row=1, col=1)
        z = np.empty((4))
        i = 0
        for yi in range(2):
            for xi in range(2):
                z[i] = fit[0] * lim_xyz[0][xi] + fit[1] * lim_xyz[1][yi] + fit[2]
                i += 1
        plane_trace = go.Mesh3d(x=[lim_xyz[0][0], lim_xyz[0][1], lim_xyz[0][0], lim_xyz[0][1]],
                                y=[lim_xyz[1][0], lim_xyz[1][0], lim_xyz[1][1], lim_xyz[1][1]],
                                z=z,
                                color = 'green',
                                opacity = 0.5,
                                delaunayaxis = 'x')
        fig.add_trace(plane_trace, row=1, col=1)
        fig.update_scenes(xaxis_range=[lim_xyz[0][0], lim_xyz[0][1]],
                          xaxis_title_text=xyz[order[0]],
                          yaxis_range=[lim_xyz[1][0], lim_xyz[1][1]],
                          yaxis_title_text=xyz[order[1]],
                          zaxis_range=[min([min(z), lim_xyz[2][0]]), max([max(z), lim_xyz[2][1]])],
                          zaxis_title_text=xyz[order[2]],
                          row=1,
                          col=1)
        fig.update_layout(title="Axis titles are good")
        fig.show()
    return fit, order


#def register(subparsers):
#    parser = subparsers.add_parser('pdbMerger',
#                                   help='Provide two .pdb files to this function to merge these into single .pdb or .gro or other supported formats.')
#    parser = addArguments(parser)
#    parser.set_defaults(func=main)


def parseArguments():
    parser = argparse.ArgumentParser(prog='pdbMerger.py',
                                     description='Provide two .pdb files to this function to merge these into single .pdb or .gro or other supported formats.',
                                     epilog='')
    parser = addArguments(parser)
    args = parser.parse_args()
    return args


def validateArguments(args):
    if len(args.segid) < 4:
        args.segid = args.segid.ljust(4)
    else:
        args.segid = args.segid[:4]
    


def addArguments(parser):
    parser.add_argument('filename1')           # positional argument
    parser.add_argument('filename2')           # positional argument
    parser.add_argument("-o",
                        "--output",
                        default="merged.pdb",
                        type=str,
                        nargs='?',
                        help='name of the output which defaults to `output.pdb`')
    parser.add_argument('-x',
                        '--x',
                        type=float,
                        nargs='?',
                        help='position along the x-axis at which center of mass should be positioned.')
    parser.add_argument('-y',
                        '--y',
                        type=float,
                        nargs='?',
                        help='position along the y-axis at which center of mass should be positioned.')
    parser.add_argument('-z',
                        '--z',
                        type=float,
                        nargs='?',
                        help='position along the z-axis at which center of mass should be positioned.')
    parser.add_argument('-box',
                        '--box',
                        type=list,
                        nargs=3,
                        help='size of the box which should be set using -box 10 10 10')
    parser.add_argument('-d',
                        '--distance',
                        type=float,
                        nargs=1,
                        help='distance from the center along the Z axis.')
    parser.add_argument('-id',
                        '--newChainID',
                        type=str,
                        default='X',
                        help='character to put in column 22 in the new pdb file corresponding to the chainID.')
    parser.add_argument('-me',
                        '--takeMe',
                        type=str,
                        help='the chainID that should be taken from the inputfile.')
    parser.add_argument('-mb',
                        '--makeBox',
                        action='store_true',  # false when not set
                        help='when set, a box will be created that is big enough to contain all atoms.')
    parser.add_argument('-zinc',
                        '--zinc',
                        type=float,
                        help='increases the z-size of the box and recenters all atoms after dealing with -box.')
    parser.add_argument('-zinc1',
                        '--zincMono',
                        action='store_true',
                        help='only increase size of 1 side of the box')
    parser.add_argument('-center',
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
    parser.add_argument('-cc',
                        '--concentric',
                        action="store_true",
                        help='position the molecule com at the main main center of mass')
    parser.add_argument('-prd',  # point residue down
                        '--thisResIDDown',
                        type=int,
                        help='specify a resID of the molecule to orient downwards the vector form molecule center of mass to residue center of mass.')
    parser.add_argument('-a',
                        '--angle',
                        type=int,
                        help='only does something if -prd/--thisResIDDown is also spcified. Obliques the plane passing throught the chain by the number of degrees specified.')
    parser.add_argument('-dbetween',
                        '--distance_between',
                        type=float,
                        help='minimum distance between minimum z value of the chain and maximum of the membrane, can be overridden by -distance.')
    parser.add_argument('-addids',
                        '--addChainIDs',
                        action="store_true",
                        help='assign chain IDs to atoms when they dont have em, as long as capital letters are not all taken. Distinct resnames get distinct chainids.')
    parser.add_argument("-sid",
                        "--segid",
                        type=str,
                        default='    ',
                        help="used to set the segment ID in column 73-76 of the .pdb output file.")
    parser.add_argument("-q",
                        "--quiet",
                        action="store_true",
                        help="suppress warnings")
    return parser


def distanceBetweenZoffset(cloud1, cloud2, target_dist=None):
    """
    Both positional arguments are n*3 (column) np.arrays
    Assumes `cloud1[2]` and `cloud2[2]` are the z-coordinates
    If target_dist is specified the function will return the clouds translated to be at this target distance from eachother
    Assumes the distance between all neighbouring points is smaller than target_dist and they are atoms at most 1 nm from eachother
    Uses only the points 10 nm down and up from the min and max z
    The same goes for offset from the cloud (as a square) borders
    Can move cloud2 towards cloud1
        cloud2 should be smaller in the xyz plane
    """

    # This is the one that has to be memoized
    @njit
    def old_squaresOfCloudsXYZ(cloud1, cloud2):
        squares = np.empty((3, cloud1.shape[0], cloud2.shape[0]))
        for d in range(3):
            for i1 in range(cloud1.shape[0]):
                for i2 in range(cloud2.shape[0]):
                    squares[d][i1][i2] = (cloud1[i1][d] - cloud2[i2][d]) ** 2
        return squares
   

    def squaresOfCloudsXYZ(cloud1, cloud2):
        cache = [{}, {}, {}]#(None, None), (None, None), (None, None)]
       
        
        @njit
        def squareDiff(col1, col2):
            returner = np.empty((col1.shape[0], col2.shape[0]))
            for i1 in range(col1.shape[0]):
                for i2 in range(col2.shape[0]):
                    returner[i1][i2] = (col1[i1] - col2[i2]) ** 2
            return returner
      
        
        def innerSquaresOfCloudsXYZ(cloud1, cloud2):
            squares = np.empty((3, cloud1.shape[0], cloud2.shape[0]))
            for d in range(3):
                serialized = np.hstack((cloud1[:,d], cloud2[:,d])).tobytes()
                if serialized in cache[d]:
                    squares[d] = cache[d][serialized]
                else:
                    squares[d] = squareDiff(cloud1[:,d], cloud2[:,d])
                    cache[d][serialized] = squares[d]
            return squares
        
        return innerSquaresOfCloudsXYZ(cloud1, cloud2)
                    
    @njit
    def distanceFromSquares(squares):
        distances = np.sqrt(np.sum(squares, axis=0))
        return np.min(distances)


    rc = 10  # Used to determine relevant range
    min_stepsize = 0.001  # used to move cloud2 closer
    dis_tol = 0.001
    step = 0
    step_tol = 1000
    negated = False
    cloud1name = "cloud1"
    cloud2name = "cloud2"
    clouddict = {cloud1name: cloud1, cloud2name: cloud2}
    centers_of_geom = {}
    mindist, maxdist = float, float
    for key in clouddict:
        centers_of_geom[key] = np.array([np.mean(clouddict[key], 0)]).flatten()
    minmax = {cloud1name: {"min": np.zeros((3)), "max": np.zeros((3))},
              cloud2name: {"min": np.zeros((3)), "max": np.zeros((3))}}
    for key, val in clouddict.items():
        for d in range(3):
            minmax[key]['min'][d] = np.min(val[:, d])
            minmax[key]['max'][d] = np.max(val[:, d])
    maxdist = minmax[cloud2name]['min'][2] - minmax[cloud1name]['max'][2]
    useful_cloud1 = []  # We dont want to calculate distances between all points
    for point in cloud1:
        if point[0] > minmax[cloud1name]['min'][0] - rc \
        and point[0] < minmax[cloud1name]['max'][0] + rc \
        and point[1] > minmax[cloud1name]['min'][1] - rc \
        and point[1] < minmax[cloud1name]['max'][1] + rc:
            if point[2] > minmax[cloud1name]['max'][2] - 10:
                useful_cloud1.append(point)
    useful_cloud1 = np.array(useful_cloud1)
    if target_dist is not None:
        distance_moved = 0
        return_cloud2 = np.array(cloud2)
        if maxdist < 0:
            return_cloud2[:,2] += target_dist - maxdist
            distance_moved += target_dist - maxdist
        mindist = distanceFromSquares(squaresOfCloudsXYZ(useful_cloud1, return_cloud2))
        while mindist - target_dist > dis_tol and step < step_tol:
            move = max([min_stepsize, mindist - target_dist])
            return_cloud2[:,2] -= move
            distance_moved     -= move
            mindist = distanceFromSquares(squaresOfCloudsXYZ(useful_cloud1, return_cloud2))
            step += 1
        print(step)
        if step == step_tol:
            print("The function distanceBetweenZoffset did not converge.", file=sys.stderr)
        return return_cloud2, distance_moved
    else:
        mindist = distanceFromSquares(squaresOfCloudsXYZ(useful_cloud1, cloud2))
        return return_cloud2, mindist


args = parseArguments()

xyz = ['x', 'y', 'z']

u_main = mda.Universe(args.filename1) # the membrane
u2 = mda.Universe(args.filename2) # the molecule


for segment in u_main.segments:
    segment.segid = args.segid

virtualbox = getVirtBox(u_main)


# Only add box dimensions when specified and big enough
if args.box is not None:
    if args.box[0] and args.box[1] and args.box[2]:
        try:
            if float(args.box[0]) >= virtualbox[0] and float(args.box[1]) >= virtualbox[1] and float(args.box[2]) >= virtualbox[2]:
                transform = boxdimensions.set_dimensions(args.box + [90, 90, 90])
                u_main.trajectory.add_transformations(transform)
            else:
                raise
        except ValueError:
            sys.exit("The specified boxdimensions are not big enought to encapsulate the system")


# Make sure the molecule in the new pdb has a different chainID
u_main_ids = getChainIDs(u_main)
chainids = getChainIDs(u2)
chainToMove = str
if args.takeMe is not None:
    if args.takeMe in chainids:
        chainToMove = args.takeMe
    else:
        print(f"The file contained no molecules with chainID {args.takeMe} so atoms in chain {chainids[0]} were taken.")
        chainToMove = chainids[0]
else:
    chainToMove = chainids[0]
u2_chain = mda.Merge(u2.select_atoms(f"chainID {chainToMove}")) # create new universe
if chainToMove in u_main_ids:
    if args.newChainID not in u_main_ids:
        for atom in u2_chain.atoms:
            atom.chainID = args.newChainID
    else:
        sys.exit(f"specify a correct chainID to assign to the molecule in {args.filename2}.")


# Add chain ids
if args.addChainIDs:  # Writer assigns if unspecified
    if not hasattr(u_main.atoms, 'chainIDs'):
        u_main.add_TopologyAttr('chainIDs')
    abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    distinct_residues_without_chainid = []
    for residue in u_main.residues:
        chainid = ""
        for atom in u_main.select_atoms("resname " + str(residue.resname)):
            if hasattr(atom, "chainID") and atom.chainID is not None and atom.chainID != " " and atom.chainID != "":
                chainid = atom.chainID
                break
        if chainid != "":
            for atom in u_main.select_atoms("resname " + str(residue.resname)):
                atom.chainID = chainid
        else:
            if residue.resname not in distinct_residues_without_chainid:
                distinct_residues_without_chainid.append(residue.resname)
    i = 0
    for nochainid in distinct_residues_without_chainid:
        while i < len(abc) and abc[i] in u_main_ids or abc[i] in chainids:
            i += 1
        if i == len(abc):
            sys.exit(f"Tried to assign more chainIDs than available in the `{abc}`")
        for atom in u_main.select_atoms("resname " + nochainid):
            atom.chainID = abc[i]
        i += 1


# Move molecule along axes relative to center of mass of u_main
com_main = u_main.atoms.center_of_geometry()
com2 = u2_chain.atoms.center_of_geometry()
if args.concentric:
    for d in range(0, 3):
        u2_chain.coord.positions[:, d] += (com_main[d] - com2[d])
        com2[d] = com_main[d]
if args.x is not None:
    offset = (args.x - com2[0] + com_main[0])  # basically new x
    u2_chain.coord.positions[:, 0] += offset
    com2[0] += offset
if args.y is not None:
    offset = (args.y - com2[1] + com_main[1])  # basically new y
    u2_chain.coord.positions[:, 1] += offset
    com2[1] += offset
if args.z is not None:
    offset = (args.z - com2[2] + com_main[2])  # basically new z
    u2_chain.coord.positions[:, 2] += offset
    com2[2] += offset


# rotate the molecule in space
ddist = 2
if args.thisResIDDown is not None:
    if args.thisResIDDown not in u2_chain.residues.resids:
        sys.exit(f"residue number {args.thisResIDDown} is not present in chain {chainToMove} in {args.filename2}.")
    res_group = u2_chain.select_atoms(f"resid {args.thisResIDDown}")
    cog_res = res_group.center_of_geometry()
    u2_chain.coord.positions -= com2
    res_group = u2_chain.select_atoms(f"resid {args.thisResIDDown}")
    cog_res = res_group.center_of_geometry()
    down_vec = np.zeros(3)
    down_vec[ddist] = -1
    angle = angleBetween(down_vec, cog_res)
    # create matrix that rotates residue towards down_vec
    if angle != 0:
        axis = np.cross(cog_res, down_vec)
        rot_mat = rotationMatrix(axis, angle)
        u2_chain.coord.positions = np.matmul(u2_chain.coord.positions, rot_mat)
    if args.angle is not None:
        if args.angle != 0:
            vecs_res_cogs = np.empty((len(u2_chain.residues), 3))  # use for cyclic peptide plane finding
            for i, residue in enumerate(u2_chain.residues):
                vecs_res_cogs[i] = u2_chain.select_atoms(f"resid {residue.resid}").center_of_geometry()
            # fit[0] * cloud[order[0]] + fit[1] * cloud[order[1]] + fit[2] ~= cloud[order[2]]
            fit, order = FindLstsqPlane(vecs_res_cogs, axis=2)  # Plane does not necessarily go through the origin
            vecs_in_plane = np.zeros((2, 3))  # use for normal_to_plane finding
            for i in range(2):
                vecs_in_plane[i][i] = 1
                vecs_in_plane[i][2] = fit[0] * vecs_in_plane[i][0] + fit[1] * vecs_in_plane[i][1]
            backorder = np.array(range(3))
            for i, d in enumerate(order):
                backorder[d] = i
            # order indexes assignment of xyz into fit
            # backorder indexes assignment of fit into xyz
            normal_to_plane = np.cross(vecs_in_plane[0], vecs_in_plane[1])
            normal_to_plane = np.array([normal_to_plane[backorder[0]], normal_to_plane[backorder[1]], normal_to_plane[backorder[2]]])
            rot_axis = np.cross(normal_to_plane, down_vec)
            res_group = u2_chain.select_atoms(f"resid {args.thisResIDDown}")
            cog_res = res_group.center_of_geometry()
            u2_chain.coord.positions -= cog_res
            ### Added because using the plane is better than just using the vector from chain cog to res cog (which will stay at the origin)
            angle = angleBetween(normal_to_plane, cog_res)
            rot_mat = rotationMatrix(rot_axis, angle - math.pi / 2)
            u2_chain.coord.positions = np.matmul(u2_chain.coord.positions, rot_mat)
            ###
            rot_mat = rotationMatrix(rot_axis, args.angle / 180 * math.pi)
            u2_chain.coord.positions = np.matmul(u2_chain.coord.positions, rot_mat)
            u2_chain.coord.positions += cog_res
    u2_chain.coord.positions += com2


# Move chain along the z-axis
distance = 0
moveIt = False
if args.distance_between is not None:
    start_time = time.time()
    u2_chain.coord.positions, distance = distanceBetweenZoffset(u_main.coord.positions, u2_chain.coord.positions, target_dist=args.distance_between)
    end_time = time.time()
    print(f"time taken to move the molecule: {end_time - start_time}")
elif args.distance is not None:
    u2_chain.coord.positions[:, ddist] += distance


# Create a new universe
u_new = mda.Merge(u2_chain.atoms, u_main.atoms)
virtualbox = getVirtBox(u_new)
newbox = []
read_box_dims = False  # Set to true when succesfully read box_dims from the input files
if args.box is not None and args.box[0] and args.box[1] and args.box[2]:
    if float(args.box[0]) >= virtualbox[0] and float(args.box[1]) >= virtualbox[1] and float(args.box[2]) >= virtualbox[2]:
        newbox = round(args.box, 3)
    else:
        sys.exit(f"Adding the molecule inside {args.filename2} caused the sepecified boxdimensions to not be big enough.")
# Only make makeshift box when set
elif args.makeBox:
    u_new.coord.positions[:, 2] -= min(u_new.coord.positions[:, 2])
    newbox = getVirtBox(u_new)
else:
    try:
        newbox = list(u_main.dimensions[0:3])
        read_box_dims = True
        # Check not all just set to 1
        num_set_to_one = 0
        for i in range(3):
            if newbox[i] - 1 < 1e-6:
                num_set_to_one += 1
        if num_set_to_one == 3:
            raise ValueError
    except (TypeError, ValueError):
        newbox = virtualbox


# modify the box?
if args.zinc is not None:
    newbox[2] += round(args.zinc, 3)
    if not args.zincMono:
        u_new.coord.positions[:, 2] += round(args.zinc / 2, 3)

# center in box
if args.center:
    com_new = [round(num, 3) for num in u_new.atoms.center_of_mass()]
    for d in range(0, 3):
        u_new.coord.positions[:, d] += (newbox[d] / 2 - com_new[d])


# make all coord positive
if args.allpos:
    for d in range(0, 3):
        mind = min(u_new.coord.positions[:, d])
        if mind < 0:
            u_new.coord.positions[:, d] -= mind

# Output
transform = boxdimensions.set_dimensions(newbox + [90, 90, 90])
u_new.trajectory.add_transformations(transform)
# Write output

context = warnings.catch_warnings() if args.quiet else contextlib.nullcontext()
with context:
    if args.quiet:
        warnings.filterwarnings("ignore", category=UserWarning, module=r"MDAnalysis\.coordinates\.PDB")
    with mda.Writer(args.output) as outfile:
        outfile.write(u_new)

if args.distance is not None or args.distance_between is not None:
    print(f"\nThe molecule in {args.filename2} was moved along the {xyz[ddist]}-axis.")

print(f"\nThe new system was stored in `{args.output}`.\n")













