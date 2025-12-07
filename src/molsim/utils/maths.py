import numpy as np
import math
from numba import njit
from scipy.linalg import lstsq


def angleBetween(v1: np.ndarray, v2: np.ndarray):
    """
    Returns the angle between 2 vectors
    """
    v1_u = v1 / np.linalg.norm(v1)  # normalizes v1
    v2_u = v2 / np.linalg.norm(v2)  # normalizes v2
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


@njit
def rotationMatrix(axis: np.ndarray, theta):
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


@njit
def cartToCyl(cloud: np.ndarray):
    """
    cloud should be a numpy.array of shape (n,3)
    """

    new_cloud = np.empty_like(cloud)
    assert len(cloud.shape) == 2
    assert cloud.shape[1] == 3
    new_cloud[:, 0] = np.sqrt(
        np.power(cloud[:, 1], 2) + np.power(cloud[:, 0], 2))
    new_cloud[:, 1] = np.arctan2(cloud[:, 1], cloud[:, 0])
    new_cloud[:, 2] = cloud[:, 2]
    return new_cloud


@njit
def cylToCart(cloud: np.ndarray):
    """
    cloud should be a numpy.array of shape (n,3)
    """
    new_cloud = np.empty_like(cloud)
    assert len(cloud.shape) == 2
    assert cloud.shape[1] == 3
    new_cloud[:, 0] = cloud[:, 0] * np.cos(cloud[:, 1])
    new_cloud[:, 1] = cloud[:, 0] * np.sin(cloud[:, 1])
    new_cloud[:, 2] = cloud[:, 2]
    return new_cloud


@njit
def distConseqPoints(cloud):
    """
    cloud should be a numpy.array of shape (n,3)

    returns a numpy.array of shape (n-1)
    """
    assert len(cloud.shape) == 2
    assert cloud.shape[1] == 3
    dist = np.empty(cloud.shape[0]-1)
    for i, vec in enumerate(cloud):
        if i == 0:
            continue
        dist[i-1] = np.linalg.norm(vec - cloud[i-1])
    return dist


def findLstsqPlane(cloud: np.ndarray, order=None, axis=None):
    """
    cloud should be a numpy.array of shape (n,3)
    order is a three element ordered permutation of [0,1,2] where in terms of dimensions fit(order[0,1]) = order[2]
    axis is only used if order is not specified
    axis tells which axis should be in the plane 0=x 1=y 2=z
    the returned order is which columns are used at the respective index
    fit[0] * cloud[order[0]] + fit[1] * cloud[order[1]] + fit[2] ~= cloud[order[2]]

    order is needed for precision ???
    """
    standard_order = [0, 1, 2]
    others = [1, 2]
    check = []
    eMes = " variable was not specified correctly."
    assert isinstance(cloud, np.ndarray) and len(
        cloud.shape) == 2 and cloud.shape[1] == 3, "The `cloud`" + eMes
    if order is not None:
        assert '__len__' in dir(order) and len(order) == 3, "The `order`" + eMes
        for d in order:
            assert d in standard_order, "The `order`" + eMes
            check.append(d)
    elif axis is not None:
        assert axis in standard_order, "The `axis`" + eMes
        order = []
        order.append(axis)
        i = 1
        for d in standard_order:
            if axis != d:
                others[i - 1] = d
                i += 1
        if max(cloud[:, others[0]]) - min(cloud[:, others[0]]) < max(cloud[:, others[1]]) - min(cloud[:, others[1]]):
            order.append(others[0])
            order.append(others[1])
        else:
            order.append(others[1])
            order.append(others[0])
    else:
        order = standard_order
    A = np.hstack((cloud[:, [order[0], order[1]]],
                  np.ones((cloud.shape[0], 1))))
    b = cloud[:, order[2]]
    fit, residual, rnk, s = lstsq(A, b)
    return fit, order

