"""
Miscellaneous helper functions.
"""

import numpy as np
import math, re
from io import BytesIO

#### python 3.6 or earlier doesn't have importlib.resources, but it's backported as importlib_resources
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

from . import data  # relative-import the *package* containing the templates
import cctk

"""
This code populates ELEMENT_DICTIONARY and ISOTOPE_DICTIONARY from a static datafile.
"""
ELEMENT_DICTIONARY = {}
ISOTOPE_DICTIONARY = {}

with pkg_resources.open_text(data, "isotopes.csv") as isotope_file:
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
    Gets atomic number from a given element symbol (converted to titlecase using ``string.title()``).

    Args:
        atomic_symbol (str): the two-character symbol

    Returns:
        the atomic number
    """
    if atomic_symbol.title() in INV_ELEMENT_DICTIONARY:
        return int(INV_ELEMENT_DICTIONARY[atomic_symbol.title()])
    else:
        raise ValueError("unknown atomic symbol: ", atomic_symbol)


"""
This code populates COVALENT_RADII_DICTIONARY from a static datafile.
"""
COVALENT_RADII_DICTIONARY = {}
with pkg_resources.open_text(data, "covalent_radii.csv") as covalent_radii:
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

"""
This code populates VDW_RADII_DICTIONARY from a static datafile.
"""
VDW_RADII_DICTIONARY = {}
with pkg_resources.open_text(data, "vdw_radii.csv") as vdw_radii:
    for line in vdw_radii:
        line_fragments = line.split(",")

        #### There's a variable number from line to line, but the first three are always number, symbol, radius
        VDW_RADII_DICTIONARY[line_fragments[0]] = line_fragments[1]

def get_vdw_radius(atomic_number):
    """
    Gets the van der Waals radius for a given element.

    Args:
        atomic_number (int): the number of the given element

    Returns:
        the van der Waals radius in Angstroms (float)
    """
    #    if isinstance(atomic_number, int):
    atomic_number = str(atomic_number)
    if atomic_number in VDW_RADII_DICTIONARY:
        return float(VDW_RADII_DICTIONARY[atomic_number])
    else:
        raise ValueError("no van der Waals radius defined for atomic number ", atomic_number)

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
        axis (np.ndarray): the vector to rotate about
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

def compute_RMSD(geometry1, geometry2, checks=True):
    """
    Computes the root mean squared difference between two geometries.

    Args:
        geometry1 (np.array (dimensions: n atoms x 3): geometry
        geometry2 (np.array (dimensions: n atoms x 3): geometry
        checks (bool): whether to check that the inputs make sense (True by default)

    Returns:
        the root-mean-square distance between the two geometries
    """
    if checks and not isinstance(geometry1, cctk.OneIndexedArray):
        raise ValueError(f"expected cctk.OneIndexedArray but got {str(type(geometry1))} instead")
    if checks and not isinstance(geometry2, cctk.OneIndexedArray):
        raise ValueError(f"expected cctk.OneIndexedArray but got {str(type(geometry2))} instead")

    if checks and len(geometry2) != len(geometry1):
        raise ValueError("can't compare two geometries with different lengths!")

    return np.sqrt( np.sum( ( geometry1.view(np.ndarray) - geometry2.view(np.ndarray) ) ** 2) / len(geometry1) )

def get_isotopic_distribution(z):
    """
    For an element with number ``z``, returns two ``np.ndarray`` objects containing that element's weights and relative abundances.

    Args:
        z (int): atomic number

    Returns:
        masses (np.ndarray): list of isotope masses
        weights (np.ndarray): list of weights (relative to 1.00 for largest)
    """
    z = str(z)
    masses = list(ISOTOPE_DICTIONARY[z].keys())
    weights = list(ISOTOPE_DICTIONARY[z].values())
    return np.array(masses), np.array(weights)

def get_avg_mass(z):
    """
    For an element with number ``z``, return average mass of that element.
    """
    masses, weights = get_isotopic_distribution(z)
    return np.dot(masses, weights)

def get_z_from_mass(desired_mass, tolerance=0.001):
    """
    For an element with atomic mass ``desired_mass``, return the element's atomic number.

    Returns ``None`` if no mass within ``tolerance`` is found.
    """
    for z in ISOTOPE_DICTIONARY.keys():
        z = int(z)
        mass = get_avg_mass(z)
        if mass == 0:
            continue

        if abs(desired_mass - mass) < tolerance:
            return z

def draw_isotopologue(z):
    """
    For an element with number ``z``, return a weighted random atomic mass (so will return 12 99% of the time and 13 1% of the time for carbon).
    """
    z = str(z)
    masses, weights = get_isotopic_distribution(z)
    return np.random.choice(masses, p=weights)

# dict: atomic symbol --> (slope, intercept)
# defines the slope to be positive
DEFAULT_NMR_SCALING_FACTORS = {
        "H" : (1.0716,  31.6660),
        "C" : (1.0300, 180.4300),
        "N" : (0.9776, 244.5626)
}

def scale_nmr_shifts(ensemble, symmetrical_atom_numbers=None, scaling_factors="default", property_name="isotropic_shielding"):
    """
    Apply linear scaling to isotropic shieldings to get chemical shifts.
    Shifts are calculated as (intercept-shielding)/slope.
    If there are no shifts available for a structure, None will be placed in both
    return lists.

    Args:
        ensemble: an ``Ensemble`` with calculated nmr shifts
        symmetrical_atom_numbers: None to perform no symmetry-averaging, a list of lists
                                  of 1-indexed atom numbers (e.g. [ [2,4,5], [7,8] ]) for
                                  a ConformationalEnsemble, or triply-nested lists for an
                                  Ensemble, where the outer index refers to the index of
                                  the Ensemble.
        scaling_factors: "default" to use DEFAULT_NMR_SCALING_FACTORS or a dict
                         (atomic symbol --> (slope,intercept)).  Elements for
                         which scaling factors are not provided will be ignored.
        property_name:   the key in properties_dict to use to locate the predicted
                         isotropic shieldings (default="isotropic_shielding")

    Returns:
        scaled_shifts: np.array (matching the shape of the original shieldings minus symmetry averaging)
        shift_labels: np.array (also matches shape)
    """
    # check inputs
    assert isinstance(ensemble, cctk.Ensemble), f"expected Ensemble but got {str(type(ensemble))} instead"
    assert len(ensemble) > 0, "empty ensemble not allowed"
    if symmetrical_atom_numbers is None:
        symmetrical_atom_numbers = []
    assert isinstance(symmetrical_atom_numbers, list), f"symmetrical atom numbers should be specified as a list of lists, but got {str(type(ensemble))} instead"
    for l in symmetrical_atom_numbers:
        assert isinstance(l, list), f"symmetrical atom numbers must be specified as lists, but got {str(type(l))} instead: {str(l)}"
    if scaling_factors == "default":
        scaling_factors = DEFAULT_NMR_SCALING_FACTORS
    else:
        assert isinstance(scaling_factors, dict)
        assert len(scaling_factors) > 0, "must provide scaling factors"
    assert isinstance(property_name, str) and len(property_name)>0, f"property_name {property_name} is invalid"

    # get shieldings and scale
    all_scaled_shifts = []
    all_shift_labels = []
    for i,(molecule,properties) in enumerate(ensemble.items()):
        if property_name in properties:
            # get atom numbers and atomic elements as OneIndexedArrays
            atomic_numbers = molecule.atomic_numbers
            n_atoms = len(atomic_numbers)
            atomic_symbols = [ get_symbol(n) for n in atomic_numbers ]
            atomic_symbols = cctk.OneIndexedArray(atomic_symbols)
            atom_numbers = list(range(1,n_atoms+1))
#            symbol_dict = dict(zip(atomic_numbers,atomic_symbols))
            all_labels = [ f"{current_symbol}{atom_number}" for current_symbol,atom_number in zip(atomic_symbols,atom_numbers) ]
            all_labels = cctk.OneIndexedArray(all_labels)

            # check symmetrical atom numbers make sense
            n_atoms = len(atomic_numbers)
            symmetrical_groups_dict = {}    # symbol --> [ [list1], [list2], ...] where each list is a group of symmetrical atom numbers
            symmetrical_groups_dict2 = {}   # symbol --> [ union of all symmetrical atom numbers for this symbol ]
#            unique_atoms_dict = {}          # symbol --> [ union of all unique atom numbers for this symbol ]
            for symmetrical_group in symmetrical_atom_numbers:
                assert len(symmetrical_group) > 1, "must be at least 2 symmetrical nuclei in a group"
                assert len(symmetrical_group) == len(set(symmetrical_group)), f"check for duplicate atom numbers in {symmetrical_group}"
                symmetrical_symbol = None
                for atom_number in symmetrical_group:
                    assert 1 <= atom_number <= n_atoms, f"atom number {atom_number} is out of range"
                    if symmetrical_symbol is None:
                        symmetrical_symbol = atomic_symbols[atom_number]
                        assert symmetrical_symbol in scaling_factors, f"no scaling factors available for the element {symmetrical_symbol}"
                    assert atomic_symbols[atom_number] == symmetrical_symbol,\
                           (f"all atoms in a symmetrical group must correspond to the same element\n"
                            f"expected element {symmetrical_symbol} for atom {atom_number},"
                            f"but got element {atomic_symbols[atom_number]}")
                if symmetrical_symbol not in symmetrical_groups_dict:
                    symmetrical_groups_dict[symmetrical_symbol] = []
                symmetrical_groups_dict[symmetrical_symbol].append(symmetrical_group)
                if symmetrical_symbol not in symmetrical_groups_dict2:
                    symmetrical_groups_dict2[symmetrical_symbol] = []
                symmetrical_groups_dict2[symmetrical_symbol].extend(symmetrical_group)

            # get shieldings
            all_shieldings = properties[property_name]

            # iterate through requested elements
            molecule_shifts = []
            molecule_labels = []
            for symbol_of_interest,(slope,intercept) in scaling_factors.items():
                # sanity checks
                assert isinstance(slope,float), f"expected slope to be float, but got {str(type(slope))}"
                assert slope != 0, "zero slope not allowed"
                assert isinstance(intercept,float), f"expected intercept to be float, but got {str(type(intercept))}"

                # determine unique atoms 
                unique_atom_numbers_list = []
                for atomic_symbol,atom_number in zip(atomic_symbols,atom_numbers):
                    if atomic_symbol != symbol_of_interest:
                        continue
                    if symbol_of_interest in symmetrical_groups_dict2:
                        if atom_number in symmetrical_groups_dict2[symbol_of_interest]:
                            continue
                    unique_atom_numbers_list.append(atom_number)

                # extract relevant shieldings and labels for unique atoms
                if len(unique_atom_numbers_list) > 0:
                    selected_shieldings = list(all_shieldings[unique_atom_numbers_list])
                    selected_labels = list(all_labels[unique_atom_numbers_list])
                else:
                    selected_shieldings = []
                    selected_labels = []

                # extract relevant shieldings and labels for symmetrical groups
                symmetrical_groups = []
                if symbol_of_interest in symmetrical_groups_dict:
                    symmetrical_groups = symmetrical_groups_dict[symbol_of_interest]
                for symmetrical_group in symmetrical_groups:
                    first_atom_number = symmetrical_group[0]
                    current_atomic_symbol = atomic_symbols[first_atom_number]
                    if current_atomic_symbol == symbol_of_interest:
                        group_shieldings = all_shieldings[symmetrical_group]
                        averaged_shielding = group_shieldings.mean()
                        selected_shieldings.append(averaged_shielding)
                        label = f"{current_atomic_symbol}"
                        for j,atom_number in enumerate(symmetrical_group):
                            label += f"{atom_number}"
                            if j < len(symmetrical_group) - 1:
                                label += "/"
                        selected_labels.append(label)

                # apply scaling
                assert len(selected_shieldings) == len(selected_labels), "shieldings and labels should have 1:1 correspondence"
                selected_shifts = np.array(selected_shieldings)
                selected_shifts = (intercept-selected_shifts)/slope
                selected_labels = np.array(selected_labels)

                # update results
                molecule_shifts.extend(selected_shifts)
                molecule_labels.extend(selected_labels)

            # update master results if appropriate
            if len(molecule_shifts) > 0:
                all_scaled_shifts.append(molecule_shifts)
                all_shift_labels.append(molecule_labels)
            else:
                # assume this means a bug
                raise ValueError("no relevant shieldings were extracted for this molecule!")
        else:
            # there are no shieldings available, so append None
            all_scaled_shifts.append(None)
            all_shift_labels.append(None)

    # return result
    scaled_shifts = np.array(all_scaled_shifts)
    shift_labels = np.array(all_shift_labels)
    return scaled_shifts, shift_labels

def compute_chirality(v1, v2, v3, v4):
    """
    Given 4 bond vectors, returns 1 or -1 based on chirality.
    For proper Cahn–Ingold–Prelog results, vectors should be passed from highest to lowest priority; however, any predictable order will give meaningful results.

    Args:
        v1 (np.ndarray): 3D bond vector
        v2 (np.ndarray): 3D bond vector
        v3 (np.ndarray): 3D bond vector
        v4 (np.ndarray): 3D bond vector

    Returns:
        value of 1 (R by CIP) or -1 (S by CIP)
    """
    assert (isinstance(v1, np.ndarray) and len(v1) == 3), "v1 needs to be a 3-element np.ndarray!"
    assert (isinstance(v2, np.ndarray) and len(v2) == 3), "v2 needs to be a 3-element np.ndarray!"
    assert (isinstance(v3, np.ndarray) and len(v3) == 3), "v3 needs to be a 3-element np.ndarray!"
    assert (isinstance(v4, np.ndarray) and len(v4) == 3), "v4 needs to be a 3-element np.ndarray!"

    e1 = np.array([1, 0, 0])
    e3 = np.array([0, 0, 1])

    # rotate v4 so that it's pointing back!
    axis1 = np.cross(v4, e1)
    theta1 = compute_angle_between(e1, v4)

    R1 = compute_rotation_matrix(axis1, theta1)
    v1 = R1 @ v1
    v2 = R1 @ v2
    v3 = R1 @ v3
    v4 = R1 @ v4

    assert 1.0 > compute_angle_between(v4, e1), "rotating v4 failed"

    # rotate v1 so that it's pointing up!
    axis2 = v4
    theta2 = compute_angle_between(np.array([0, v1[1], v1[2]]), e3) # projection of v1 onto e2•e3 plane

    R2 = compute_rotation_matrix(axis2, np.sign(v1[1]) * theta2)
    v1 = R2 @ v1
    v2 = R2 @ v2
    v3 = R2 @ v3
    v4 = R2 @ v4

    assert 1.0 > compute_angle_between(v4, e1), f"rotating v4 failed - 1.0 ≤ {compute_angle_between(v4, e1)}"
    assert 1.0 > compute_angle_between(e3, np.array([0, 0, v1[2]])), f"rotating v1 failed, - 1.0 ≤ {compute_angle_between(e3, np.array([0, 0, v1[2]]))}"

    answer = np.sign(v2[1])
    assert np.sign(v3[1]) != answer, "at this point our two substituents are on the same side of the plane that's supposed to divide them"
    return answer

# constants for calculating entropy
ENTROPY_FACTOR_1 = 1.43877695998381562 # 2.99792458E10 * 6.62606957E-34 / 1.3806488E-23
ENTROPY_FACTOR_2 = 1.9872041348        # 8.3144621 / 4.184

def get_entropy(frequencies, temperature):
    """
        Computes the total entropy of a given set of frequencies.

        Args:
            frequencies (list): in cm-1
            temperature (float): in K

        Returns:
            entropy (float): in hartree
    """
    factor0 = ENTROPY_FACTOR_1 / temperature
    entropy = 0.0
    for frequency in frequencies:
        factor = factor0 * frequency
        temp = factor * 1.0/(math.exp(factor)-1.0) - math.log(1.0-math.exp(-factor))
        temp = temp * ENTROPY_FACTOR_2
        entropy += temp
    return entropy / 627.509469

def get_corrected_free_energy(free_energy, frequencies, frequency_cutoff=100.0, temperature=298.15):
    """
        Computes the free energy by moving all positive frequencies below ``frequency_cutoff``
        to the cutoff.  See Cramer/Truhlar, J. Phys. Chem. B, 2011, 115, 14556.

        Args:
            free_energy (float): in hartree
            frequencies (list): in cm-1
            frequency_cutoff (float): in cm-1
            temperature (float): in K

        Returns:
            corrected_free_energy (float): in hartree
    """
    low_frequencies = []
    for frequency in frequencies:
        if frequency > 0 and frequency < frequency_cutoff:
            low_frequencies.append(frequency)
    entropy_uncorrected = get_entropy(low_frequencies, temperature)
    entropy_corrected = get_entropy([frequency_cutoff], temperature) * len(low_frequencies)
    entropy_correction = (entropy_uncorrected - entropy_corrected)*temperature/1000.0
    corrected_free_energy = free_energy + entropy_correction
    return corrected_free_energy

def numpy_to_bytes(arr):
    """ Utility function for pickling numpy arrays """
    arr_bytes = BytesIO()
    np.save(arr_bytes, arr, allow_pickle=True)
    arr_bytes = arr_bytes.getvalue()
    return arr_bytes

def bytes_to_numpy(arr_bytes):
    """ Utility function for unpickling numpy arrays """
    load_bytes = BytesIO(arr_bytes)
    loaded_np = np.load(load_bytes, allow_pickle=True)
    return loaded_np

def compute_mass_spectrum(formula_dict, **kwargs):
    """
    Computes the expected low-res mass spec ions for a given formula.

    Args:
        formula dict (dict): e.g. {"C": 6, "H": 6}

    Returns:
        list of m/z ions
        list of relative weights (out of 1 total)
    """
    form_vec = np.zeros(shape=92, dtype=np.int8)
    for z, n in formula_dict.items():
        if isinstance(z, str):
            z = get_number(z)
        assert isinstance(z, int), "atomic number must be integer"
        form_vec[z] += n

    masses, weights = _recurse_through_formula(form_vec, [0], [1], **kwargs)

    new_masses, indices = np.unique(np.round(masses, decimals=1), return_inverse=True)
    new_weights = np.zeros_like(new_masses)
    for k in range(len(new_weights)):
        new_weights[k] = np.sum(weights[np.nonzero(indices == k)])
    new_weights = new_weights / np.max(new_weights)

    return new_masses, new_weights

def _recurse_through_formula(formula, masses, weights, cutoff=0.0000001, mass_precision=4, weight_precision=8):
    """
    Recurses through a formula and generates m/z isotopic pattern using tail recursion.

    To prevent blowup of memory, fragments with very low abundance are ignored. Masses and weights are also rounded after every step.
    To prevent error accumulation, internal precisions several orders of magnitude lower than the precision of interest should be employed.
    The default values should work nicely for low-res MS applications.

    Args:
        formula (np.ndarray, dtype=np.int8): vector containing atoms left to incorporate. first element should always be 0 as there is no element 0.
        masses (np.ndarray): list of mass fragments at current iteration
        weights (np.ndarray): relative weights at current iteration
        cutoff (float): cutoff for similarity (masses within ``cutoff`` will be combined)
        mass_precision (int): number of decimal places to store for mass
        weight_precision (int): number of decimal places to store for weight

    Returns:
        masses
        weights
    """
    # check how many elements we haven't recursed thru yet
    if np.array_equal(formula, np.zeros(shape=92, dtype=np.int8)):
        return masses[np.argsort(masses)], weights[np.argsort(masses)]

    # get masses/weights for current element
    current_e = np.nonzero(formula)[0][0]
    e_masses, e_weights = get_isotopic_distribution(current_e)

    # combinatorially add the new masses and weights to our current lists
    new_masses = np.zeros(shape=(len(masses)*len(e_masses)))
    new_weights = np.zeros(shape=(len(masses)*len(e_masses)))
    for i in range(len(masses)):
        for j in range(len(e_masses)):
            new_masses[i*len(e_masses)+j] = masses[i] + e_masses[j]
            new_weights[i*len(e_masses)+j] = weights[i] * e_weights[j]

    # delete duplicates and adjust weights (complicated)
    newer_masses, indices = np.unique(np.round(new_masses, decimals=mass_precision), return_inverse=True)
    newer_weights = np.zeros_like(newer_masses)
    for k in range(len(newer_weights)):
        newer_weights[k] = np.sum(new_weights[np.nonzero(indices == k)])
    newer_weights = np.round(newer_weights, decimals=weight_precision)

    # prune the low-abundance masses/weights and move on to the next element
    formula[current_e] += -1
    above_cutoff = np.nonzero(newer_weights > cutoff)
    return _recurse_through_formula(formula, newer_masses[above_cutoff], newer_weights[above_cutoff], cutoff, mass_precision, weight_precision)

def formula_dict_from_string(formula_string):
    """
    Eugene challenged me to code golf, this isn't my fault.

    Args:
        formula_string (str): the formula as a string, e.g. C10H12N2O1. you need the "1" explicitly

    Returns:
        formula_dict (dict): e.g. {'C': 10, 'H': 12, 'N': 2, 'O': 1}
    """
    return {t[0]: int(t[1]) for t in re.findall(r"([a-z]+)([0-9]+)", formula_string, re.I)}
