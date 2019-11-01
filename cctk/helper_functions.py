import numpy as np
import math
import re

#### python 3.6 or earlier doesn't have importlib.resources, but it's backported as importlib_resources
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources 

from . import data # relative-import the *package* containing the templates

"""
This code populates ELEMENT_DICTIONARY from a static datafile.
"""
ELEMENT_DICTIONARY = {}
#with open('data/isotopes.csv', mode='r') as isotope_file:
isotope_file = pkg_resources.open_text(data, 'isotopes.csv')
for line in isotope_file: 
    symbol, number, mass, abundance  = line.split(',')
    if symbol == "Symbol":
        continue
    ELEMENT_DICTIONARY[number] = symbol

INV_ELEMENT_DICTIONARY = {v: k for k, v in ELEMENT_DICTIONARY.items()}

def get_symbol(atomic_number):
    """ 
    Gets element symbol from a given atomic number.
    
    Args:
        atomic_number (int): the number of the given element

    Returns: 
        the two-character atomic symbol string
    """
    if isinstance(atomic_number, int):
        atomic_number = str(atomic_number)
    if atomic_number in ELEMENT_DICTIONARY:
        return ELEMENT_DICTIONARY[atomic_number]
    else:
        raise ValueError("unknown atomic number: ", atomic_number)

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
covalent_radii = pkg_resources.open_text(data, 'covalent_radii.csv')
for line in covalent_radii:
    line_fragments = line.split(',')

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
    
    if isinstance(atomic_number, int):
        atomic_number = str(atomic_number)
    if atomic_number in COVALENT_RADII_DICTIONARY:
        return float(COVALENT_RADII_DICTIONARY[atomic_number])
    else:
        raise ValueError("no covalent radius defined for atomic number ", atomic_number) 

def compute_distance_between(v1, v2):
    """ 
    Computes the L2 distance between two vectors.
    """
    return np.linalg.norm(v1-v2)

# normalizes the given vector so that it has unit length
def compute_unit_vector(vector):
    """
    Normalizes a vector, returning a unit vector pointing in the same direction.
    Returns the zero vector if the zero vector is given. 
    """
    if np.linalg.norm(vector) == 0:
        return vector
    else:
        return vector / np.linalg.norm(vector)

# compute the angle between two vectors in degrees
def compute_angle_between(v1, v2, unit='degree'):
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
    if unit == 'degree':
        return to_degrees(angle)
    elif unit == 'radian':
        return angle
    else: 
        raise ValueError(f"invalid unit {unit}: must be 'degree' or 'radian'!")

# compute the dihedral angle in degrees
def compute_dihedral_between(p0,p1,p2,p3):
    """
    Computes the dihedral angle between four points.
    """
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 = compute_unit_vector(b1)

    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def search_for_block(lines, start, end, count=1):
    """
    Search through a file (lines) and locate a block starting with "start" and ending with "end".
    
    Args: 
        lines (list): a list of the lines in the file
        start (str): a pattern that matches the start of the block (can contain special characters)
        end (str): a pattern that matches the end of the block (can contain special characters)
        count (int): how many matches to search for
    
    Returns: 
        a single match (str) if count == 1 or a list of matches (str) if count > 1.
    """
    current_match = '' 
    match = [None] * count

    start_pattern = re.compile(start)
    end_pattern = re.compile(end)

    index = 0
    for line in lines:
        if current_match:
            if end_pattern.search(line): 
                match[index] = current_match
                current_match = None
                index += 1

                if index == count:
                    break
            else:
                current_match = current_match + line 
        else:
            if start_pattern.search(line):
                current_match = line

    if count == 1:
        return match[0]
    else:
        return match

def compute_rotation_matrix(axis, theta):
    """
    Return the rotation matrix for rotation around `axis` by `theta` degrees..
    Adapted from user "unutbu" on StackExchange. 
    
    Args:
        axis (vector): the vector to rotate about
        theta (float): how much to rotate
    
    Returns: 
        the 3x3 rotation matrix
    """

    if (not isinstance(axis, np.ndarray)) or (len(axis) != 3):
        raise TypeError("axis must be np array with 3 elements")

    try:
        theta = float(theta)
    except:
        raise TypeError("theta must be float!")

    # convert to radians
    theta = to_radians(theta)
    axis = compute_unit_vector(axis) 
   
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def to_radians(theta):
    return (theta * math.pi) / 180

def to_degrees(theta):
    return (theta * 180) / math.pi
