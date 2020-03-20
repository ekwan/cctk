import numpy as np
import math
import re

#### python 3.6 or earlier doesn't have importlib.resources, but it's backported as importlib_resources
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from . import data  # relative-import the *package* containing the templates

"""
This code populates ELEMENT_DICTIONARY and ISOTOPE_DICTIONARY from a static datafile.
"""
ELEMENT_DICTIONARY = {}
ISOTOPE_DICTIONARY = {}
isotope_file = pkg_resources.open_text(data, "isotopes.csv")
prev_number = 1
current_dict = {}
for line in isotope_file:
    symbol, number, mass, abundance = line.split(",")
    if symbol == "Symbol":
        continue

    ELEMENT_DICTIONARY[number] = symbol

    if number == prev_number:
        current_dict[float(mass)] = float(abundance.rstrip())
    else:
        ISOTOPE_DICTIONARY[prev_number] = current_dict
        current_dict = {}
        current_dict[float(mass)] = float(abundance.rstrip())

    prev_number = number

ISOTOPE_DICTIONARY[prev_number] = current_dict
ELEMENT_DICTIONARY["0"] = "Bq"

INV_ELEMENT_DICTIONARY = {v: int(k) for k, v in ELEMENT_DICTIONARY.items()}


def get_symbol(atomic_number):
    """
    Gets element symbol from a given atomic number.

    Args:
        atomic_number (int): the number of the given element

    Returns:
        the two-character atomic symbol string
    """
    atomic_number = str(atomic_number)
    if atomic_number in ELEMENT_DICTIONARY:
        return ELEMENT_DICTIONARY[atomic_number]
    else:
        raise ValueError(f"unknown atomic number: '{atomic_number}'")


def get_number(atomic_symbol):
    """
    Gets atomic number from a given element symbol.

    Args:
        atomic_symbol (str): the two-character symbol

    Returns:
        the atomic number
    """
    if atomic_symbol in INV_ELEMENT_DICTIONARY:
        return int(INV_ELEMENT_DICTIONARY[atomic_symbol])
    else:
        raise ValueError("unknown atomic symbol: ", atomic_symbol)


"""
This code populates COVALENT_RADII_DICTIONARY from a static datafile.
"""
COVALENT_RADII_DICTIONARY = {}
covalent_radii = pkg_resources.open_text(data, "covalent_radii.csv")
for line in covalent_radii:
    line_fragments = line.split(",")

    #### There's a variable number from line to line, but the first three are always number, symbol, radius
    if line_fragments[1] == "Symbol":
        continue
    COVALENT_RADII_DICTIONARY[line_fragments[0]] = line_fragments[2]


def get_covalent_radius(atomic_number):
    """
    Gets the covalent radius for a given element.

    Args:
        atomic_number (int): the number of the given element

    Returns:
        the covalent radius in Angstroms (float)
    """
    #    if isinstance(atomic_number, int):
    atomic_number = str(atomic_number)
    if atomic_number in COVALENT_RADII_DICTIONARY:
        return float(COVALENT_RADII_DICTIONARY[atomic_number])
    else:
        raise ValueError("no covalent radius defined for atomic number ", atomic_number)


def compute_distance_between(v1, v2, _norm=np.linalg.norm):
    """
    Computes the L2 distance between two vectors.

    (preloading ``_norm`` speeds repeated calls, since Python doesn't have to look up the function every time)
    """
    return _norm(v1 - v2)


def compute_unit_vector(vector):
    """
    Normalizes a vector, returning a unit vector pointing in the same direction.
    Returns the zero vector if the zero vector is given.
    """
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector
    else:
        return vector / norm


def compute_angle_between(v1, v2, unit="degree"):
    """
    Computes the angle between two vectors.

    Args:
        v1 (ndarray): first vector
        v2 (ndarray): second vector
        unit (str): 'degree' or 'radian'

    Returns:
        the angle between the two vectors
    """
    v1_u = compute_unit_vector(v1)
    v2_u = compute_unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if unit == "degree":
        return np.degrees(angle) % 360
    elif unit == "radian":
        return angle % (2 * math.pi)
    else:
        raise ValueError(f"invalid unit {unit}: must be 'degree' or 'radian'!")


def compute_dihedral_between(p0, p1, p2, p3, unit="degree"):
    """
    Computes the dihedral angle between four points.
    """
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    b1 = compute_unit_vector(b1)

    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    angle = np.arctan2(y, x)

    if unit == "degree":
        return np.degrees(angle) % 360
    elif unit == "radian":
        return angle % (2 * math.pi)
    else:
        raise ValueError(f"invalid unit {unit}: must be 'degree' or 'radian'!")


def compute_rotation_matrix(axis, theta):
    """
    Return the rotation matrix for rotation around ``axis`` by ``theta`` degrees..
    Adapted from user "unutbu" on StackExchange.

    Args:
        axis (vector): the vector to rotate about
        theta (float): how much to rotate (in degrees)

    Returns:
        the 3x3 rotation matrix
    """
    if (not isinstance(axis, np.ndarray)) or (len(axis) != 3):
        raise TypeError("axis must be np array with 3 elements")

    try:
        theta = float(theta)
    except:
        raise TypeError("theta must be float!")

    theta = np.radians(theta)
    axis = compute_unit_vector(axis)

    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)

    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )

def align_matrices(P_partial, P_full, Q_partial, return_matrix=False):
    """
    Rotates one set of points onto another using the Kabsch algorithm.
    The rotation that best aligns P_partial into Q_partial will be found and then applied to P_full.

    Args:
        P_partial (matrix): atoms of P that correspond to Q
        P_full (matrix): full matrix to rotate
        Q (matrix): matrix to align to

    Returns:
        rotated P matrix
    """
    assert np.shape(P_partial) == np.shape(Q_partial)

    C = P_partial.T @ Q_partial
    U, S, Vt = np.linalg.svd(C)

    V = Vt.T
    d = np.linalg.det(V @ U.T)
    middle = np.identity(3)

    if d < 0.0:
        middle[2][2] = -1.0

    rotation = U @ middle @ Vt
    return P_full @ rotation


def compute_RMSD(geometry1, geometry2):
    """
    Computes the root mean squared difference between two geometries.

    Args:
        geometry1 (list, or equivalent): first geometry
        geometry2 (list, or equivalent): second geometry
    """
    if len(geometry1) != len(geometry2):
        raise ValueError("can't compare two geometries with different lengths!")

    try:
        geometry1 = np.array(geometry1)
        geometry2 = np.array(geometry2)
    except:
        raise TypeError("geometries cannot be cast to numpy arrays!")

    squared_difference = np.square(geometry1 - geometry2)
    temp = np.sum(squared_difference) / (3 * len(geometry1))
    return np.sqrt(temp)

def get_isotopic_distribution(z):
    """
    For an element with number ``z``, returns two ``np.array`` objects containing that element's weights and relative abundances.
    """
    z = str(z)
    masses = list(ISOTOPE_DICTIONARY[z].keys())
    weights = list(ISOTOPE_DICTIONARY[z].values())
    return np.array(masses), np.array(weights)
