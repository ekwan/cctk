import numpy as np
import re

ELEMENT_DICTIONARY = {}
with open('data/isotopes.csv', mode='r') as isotope_file:
    for line in isotope_file: 
        symbol, number, mass, abundance  = line.split(',')
        if symbol == "Symbol":
            continue
        ELEMENT_DICTIONARY[number] = symbol

def get_symbol(atomic_number):
    ''' 
    Get element symbol from a given atomic number.
    ''' 
    if isinstance(atomic_number, int):
        atomic_number = str(atomic_number)
    if atomic_number in ELEMENT_DICTIONARY:
        return ELEMENT_DICTIONARY[atomic_number]
    else:
        raise ValueError("unknown atomic number: ", atomic_number)

# compute the L2 distance between v1 and v2
def compute_distance_between(v1, v2):
     return np.linalg.norm(v1-v2)

# normalizes the given vector so that it has unit length
def compute_unit_vector(vector):
    return vector / np.linalg.norm(vector)

# compute the angle between two vectors in degrees
def compute_angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

# compute the dihedral angle in degrees
def compute_dihedral_between(p0,p1,p2,p3):
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

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
    '''
    Search through a file (lines) and locate a block starting with "start" and ending with "end".
    Returns a string if count == 1 or an array of strings if count > 1.
    '''
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
