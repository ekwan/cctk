import numpy as np
from glob import glob
import re

### Helper Functions ###

# element dictionary
ELEMENT_DICTIONARY = {"1":"H", "6":"C", "7":"N", "8":"O", "9":"F"}
def get_symbol(atomic_number):
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

# represents a collection of conformers or geometries
class Geometries(object):
    # geometries: np.array (conformer, x,y,z)
    # all_symbols_list: list of lists (conformer, atom symbol)
    # is_conformers: True if the geometries belong to conformers
    # bonds: list of lists of 3-tuples (1-indexed atom numbers, bond order)
    #        outer list is geometries
    # energies: parallel list of energies or None if no energies
    def __init__(self, geometries, all_symbols_list, bonds=None, energies=None):
        # check that there is at least one geometry
        if len(geometries) == 0:
            raise ValueError("empty geometries!")
        self.geometries = np.array(geometries)
        self.n_geometries = len(geometries)

        # check that there is one symbol list for every geometry
        if len(geometries) != len(all_symbols_list):
            raise ValueError("mismatch in list sizes between geometries and symbols!")
        self.all_symbols_list = all_symbols_list
        
        # store bonds
        n = len(geometries[0])
        if bonds is not None:
            for this_bonds,geometry in zip(bonds, geometries):
                for i,j,k in this_bonds:
                    if i < 1 or i > n or j < 1 or j > n or i == j or k < 1:
                        raise ValueError("invalid bond: %d-%d" % (i,j))
            self.bonds = bonds
        else:
            self.bonds = None
        
        # determine if this is a group of conformers
        self._is_conformers = True
        initial_symbol_list = all_symbols_list[0]
        for i, other_symbol_list in enumerate(all_symbols_list):
            # check that there is one symbol for every atom vector
            if len(other_symbol_list) != len(geometries[i]):
                raise ValueError("unexpected length of symbols list")
            # determine whether all symbol lists are the same
            elif initial_symbol_list != other_symbol_list:
                self._is_conformers = False
                break
        
        # store energies
        if energies is None:
            self.energies = None
        else:
            if not isinstance(energies,list) and not isinstance(energies,np.ndarray):
                raise ValueError("energies must be a list")
            if len(energies) != len(geometries):
                raise ValueError("mismatch in number of energies (%d) vs. geometries (%d)" % (len(energies), len(geometries)))
            for e in energies:
                if e is None:
                    raise ValueError("cannot have None as an energy type")
                elif not isinstance(e, float):
                    raise ValueError("got unexpected energy type")
            self.energies = np.array(energies)

    def __str__(self):
        n_geometries = len(self.geometries)
        min_atoms = len(self.geometries[0])
        if self._is_conformers:
            return "%d conformers (%d atoms)" % (n_geometries, min_atoms)
        else:
            max_atoms = len(self.geometries[0])
            for geometry in self.geometries:
                n_atoms = len(geometry)
                if n_atoms > max_atoms:
                    max_atoms = n_atoms
                elif n_atoms < min_atoms:
                    min_atoms = n_atoms
            return "%d geometries (%d-%d atoms)" % (n_geometries, min_atoms, max_atoms)
   
    # write all the geometries to a .mol2
    # filename should include the .mol2 extension
    def write_mol2(self, filename, molecule_name="untitled"):
        with open(filename, "w") as mol2:
            for i in range(self.n_geometries):
                geometry = self.geometries[i]
                symbols = self.all_symbols_list[i]
                bonds = []
                if self.bonds is not None:
                    bonds = self.bonds[i]
                n_atoms = len(geometry)
                n_bonds = len(bonds)
                print("@<TRIPOS>MOLECULE", file=mol2)
                print(molecule_name, file=mol2)
                print("%d %d" % (n_atoms, n_bonds), file=mol2)
                print("SMALL\nNO CHARGES\n", file=mol2)
                print("@<TRIPOS>ATOM", file=mol2)
                i=0
                for symbol,xyz in zip(symbols,geometry):
                    i += 1
                    x,y,z = xyz
                    print("%d %s%d %12.8f %12.8f %12.8f %s" %
                          (i, symbol, i, x, y, z, symbol), file=mol2)
                print("@<TRIPOS>BOND",file=mol2)
                i=0
                for atom1,atom2,bond_order in bonds:
                    i += 1
                    print("%d %d %d %d" % (i, atom1, atom2, bond_order), file=mol2)
                print("", file=mol2)
    
    # write the specified geometry (0-indexed) to a .gjf
    # filename should include the .gjf extension
    def write_gjf(self, geometry_index, filename, header="#\n\ntitle\n\n0 1", footer="\n\n"):
        # ensure geometry_number is valid
        if geometry_index < 0 or geometry_index >= len(self.geometries):
            raise ValueError("invalid geometry index for gjf dump")
        
        # retrieve specified geometry and atom symbols
        this_geometry = self.geometries[geometry_index]
        this_symbols = self.all_symbols_list[geometry_index]
        
        # write result to file
        with open(filename, "w") as gjf:
            print(header, file=gjf)
            for symbol, coordinates in zip(this_symbols, this_geometry):
                x,y,z = coordinates
                print("%2s %12.8f %12.8f %12.8f" % (symbol, x, y, z), file=gjf)
            print(footer, file=gjf)
            print("\n", file=gjf)

    # compute distance for geometry i (0-indexed)
    # between atoms a1 and a2 (1-indexed)
    def distance(self,i,a1,a2):
        return compute_distance_between(self.geometries[i,a1-1], self.geometries[i,a2-1])
    
    # compute angle in degrees for geometry i (0-indexed)
    # between atoms a1 and a2 and a3 (1-indexed)
    def angle(self,i,a1,a2,a3):
        v1 = self.geometries[i,a1-1] - self.geometries[i,a2-1]
        v2 = self.geometries[i,a3-1] - self.geometries[i,a2-1]
        return np.degrees(compute_angle_between(v1,v2))
    
    # compute dihedral angle in degrees for geometry i (0-indexed)
    # between atoms a1 and a2 and a3 and a4 (1-indexed)
    def dihedral(self,i,a1,a2,a3,a4):
        return compute_dihedral_between(self.geometries[i,a1-1],
                                        self.geometries[i,a2-1],
                                        self.geometries[i,a3-1],
                                        self.geometries[i,a4-1])

    # rotates one set of points onto another using Kabsch algorithm
    # the rotation that best aligns P_partial onto Q_partial will be found
    # and then applied to P_full (partial means only comparison atoms are kept)
    @staticmethod
    def _align(P_partial,P_full,Q_partial):
        assert(np.shape(P_partial) == np.shape(Q_partial))
        C = P_partial.T @ Q_partial
        U,S,Vt = np.linalg.svd(C)
        V = Vt.T
        d = np.linalg.det(V @ U.T)
        middle = np.identity(3)
        if d < 0.0:
            middle[2][2]=-1.0
        rotation = U @ middle @ Vt
        return P_full @ rotation
   
    # returns a new Geometries object where all geometries have been
    # aligned to the specified geometry number (1-indexed)
    # atom_numbers: which atoms to perform the alignment to
    # heavy_atoms: pick the non-hydrogen atoms automatically
    def align_all_geometries(self, geometry_number=1, atom_numbers = None, heavy_atoms = True):
        if not self._is_conformers:
            raise ValueError("cannot perform alignment on non-conformers")
        
        # extract and center geometries for alignment
        atom_indices = self._get_atom_indices(atom_numbers, heavy_atoms)
        if len(atom_indices) < 3:
            raise ValueError("not enough atoms for alignment")
        partial_geometries = self.geometries[:,atom_indices,:]
        partial_geometry_centroids = [ geometry.mean(axis=0) for geometry in partial_geometries ]
        partial_geometries = [ geometry - centroid for geometry,centroid in 
                               zip(partial_geometries, partial_geometry_centroids) ]
        partial_template_geometry = partial_geometries[geometry_number-1]
        full_geometries = [ geometry - centroid for geometry,centroid in
                            zip(self.geometries, partial_geometry_centroids) ]
        
        # perform alignment
        print("Aligning %d structures using %d atoms..." %
              (len(self.geometries),len(atom_indices)), end="")
        aligned_geometries = [ Geometries._align(P_partial,P_full,partial_template_geometry) for
                               P_partial,P_full in
                               zip(partial_geometries,full_geometries) ]
        print("done!")
        
        # return result
        aligned_geometries = np.array(aligned_geometries)
        new_bonds = None
        if self.bonds is not None:
            new_bonds = self.bonds.copy()
        new_energies = None
        if self.energies is not None:
            new_energies = self.energies.copy()
        return Geometries(aligned_geometries, self.all_symbols_list.copy(), new_bonds, new_energies)

    # compute RMSD between two geometries (0-indexed)
    # geometries are taken as is
    @staticmethod
    def _compute_RMSD(geometry1, geometry2):
        squared_difference = np.square(geometry1 - geometry2)
        temp = np.sum(squared_difference) / (3*len(geometry1))
        return np.sqrt(temp)

    # extract partial geometry based on a list
    def _extract_partial_geometries(self, atom_indices):
        if not self._is_conformers:
            raise ValueError("can't extract partial geometries for non-conformers")
        return self.geometries[:,atom_indices,:]
    
    # extract some of the geometries to a new object
    # geometry_indices: list of 0-indexed geometry numbers
    def extract_geometries(self, geometry_indices):
        new_geometries = self.geometries[geometry_indices].copy()
        new_symbols = [ self.all_symbols_list[i] for i in geometry_indices ]
        new_bonds = None
        if self.bonds is not None:
            new_bonds = [ self.bonds[i] for i in geometry_indices ]
        new_energies = None
        if self.energies is not None:
            new_energies = [ self.energies[i] for i in energies ]
        return Geometries(new_geometries,new_symbols,new_bonds,new_energies)
    
    # get a 0-indexed list of atom indices corresponding to the specified atoms
    # atom_numbers: an explicit of 1-indexed atoms, or None for all atoms
    # heavy_atoms:  if atom_number_list is None and this is True, then automatically
    #               populate with non-hydrogen atom indices
    def _get_atom_indices(self, atom_numbers, heavy_atoms):
        # ensure that it makes sense to get atom indices
        if not self._is_conformers:
            raise ValueError("cannot extract atom indices on non-conformers")
        
        # get atom indices
        n_atoms = len(self.geometries[0])
        atom_indices = []
        
        if atom_numbers is None:
            # explicit atom numbers not provided
            if heavy_atoms:
                # determine the indices of the heavy atoms
                for i,symbol in enumerate(self.all_symbols_list[0]):
                    if symbol != "1" and symbol != "H":
                        atom_indices.append(i)
            else:
                # provide all atom indices
                atom_indices = list(range(n_atoms))
                
        else:
            # explicit comparison atoms provided,
            # so warn if heavy_atoms is set
            if heavy_atoms:
                print("warning: explicit comparison atoms provided, so ignoring heavy_atoms flag")
            
            # check comparison atoms are valid
            for i in atom_numbers:
                # check that comparison atoms are valid
                if i < 1 or i > n_atoms:
                    raise ValueError("invalid comparison atom")
                    
            # convert atom numbers to atom indices
            atom_indices = [ i-1 for i in atom_numbers ]
            
        # return result
        return atom_indices
          
    # returns non-redundant conformations
    # when redundancies are found, only the first geometry is kept
    # assumes geometries have already been aligned
    # atom_numbers is 1-indexed; if present, overrides heavy_atoms
    def eliminate_redundant_conformers(self, RMSD_threshold, atom_numbers = None,
                                       heavy_atoms = True, print_progress=True, update_interval=500):
        # check if conformer elimination makes sense
        if not self._is_conformers:
            raise ValueError("cannot perform conformer elimination on non-conformers")
        
        # determine which atoms to compare for RMSD calculations
        comparison_atom_indices = self._get_atom_indices(atom_numbers, heavy_atoms)
        if print_progress:
            print("Eliminating redundant conformers.  Will use these atom numbers for comparison:")
            print("  ", end="")
            print_counter = 0
            for i in comparison_atom_indices:
                symbol = self.all_symbols_list[0][i]
                atom_number = i+1
                print("%s%d " % (symbol, atom_number), end="")
                print_counter += 1
                if print_counter > 15:
                    print_counter = 0
                    print("")
                    print("  ", end="")
            print()
        
        # extract partial geometries
        old_partial_geometries = self._extract_partial_geometries(comparison_atom_indices)
        
        # initialize two parallel lists
        new_partial_geometries = []          # contains unique partial geometries
        new_partial_geometry_indices = []    # corresponding indices from original list of geometries
        
        # add the first geometry to the lists
        new_partial_geometries.append( old_partial_geometries[0] )
        new_partial_geometry_indices.append(0)
        
        # initialize some counters to keep track of progress
        current_comparisons = 0
        last_comparisons = 0
        n_geometries = len(old_partial_geometries)
        total_comparisons = n_geometries*(n_geometries-1)/2
        
        # remove redundant conformers
        for counter, candidate_geometry in enumerate(old_partial_geometries[1:]):
            okay_to_add = True
            for established_geometry in new_partial_geometries:
                RMSD = Geometries._compute_RMSD(candidate_geometry, established_geometry)
                #print(counter,RMSD)
                if RMSD < RMSD_threshold:
                    okay_to_add = False
                    break
            if okay_to_add:
                #print("added")
                new_partial_geometries.append(candidate_geometry)
                new_partial_geometry_indices.append(counter+1)
            #else:
            #    print("not added")
            #print()
                    
            current_comparisons += n_geometries - counter - 1
            if print_progress and current_comparisons - last_comparisons >= update_interval:
                last_comparisons = current_comparisons
                percentage = 100.0*(current_comparisons/total_comparisons)
                print("loop %d of %d; %d of %d comparisons (%.1f%%)" %
                      (counter+1, n_geometries-1, current_comparisons, total_comparisons, percentage),
                       end="\r", flush=True)
        
        # provide final status update
        if print_progress:
            print("  loop %d of %d; %d of %d comparisons (%.1f%%)" %
                  (n_geometries-1, n_geometries-1, total_comparisons, total_comparisons, 100.0),
                  flush=True)
            print("Conformer elimination complete.  %d of %d geometries kept." %
                  (len(new_partial_geometries), len(old_partial_geometries)))
            
        # figure out which original geometries correspond to
        # the retained partial geometries and return the result
        new_complete_geometries = self.geometries[new_partial_geometry_indices,:,:]
        new_symbols_list = [ self.all_symbols_list[0].copy() for i in range(len(new_complete_geometries)) ]
        new_bonds = None
        if self.bonds is not None:
            new_bonds = [ self.bonds[0].copy() for i in range(len(new_complete_geometries)) ]
        new_energies = None
        if self.energies is not None:
            new_energies = [ self.energies[i] for i in new_partial_geometry_indices ]
        return Geometries(new_complete_geometries, new_symbols_list, new_bonds, new_energies)

### Mol2 File Format Reader ###
# reads a mol2 file and
# returns a numpy array:
# conformer, x, y, z
def read_conformers_from_mol2(filename):
    # read file
    print("Reading '%s'..." % filename, end='')
    with open(filename, 'r') as filehandle:
        lines = filehandle.read().splitlines()
    print("read %d lines..." % len(lines), end='')
    
    # initialize arrays
    all_geometries = []
    all_symbol_lists = []
    all_bonds = []
    this_geometry = []
    this_symbol_list = []
    this_bonds = []
    
    # parse file
    i=0
    in_geometry_block = False
    in_bond_block = False
    bond_number = 0
    while i < len(lines):
        # read the current line
        line = lines[i]
        
        # determine if we are in a geometry block
        if line.startswith("@<TRIPOS>ATOM"):
            # step forward to the first geometry line
            in_geometry_block = True
            in_bond_block = False
            i += 1
            line = lines[i]
        elif line.startswith("@<TRIPOS>BOND"):
            # reached the end of a geometry block, so store geometry
            in_geometry_block = False
            in_bond_block = True
            bond_number = 0
            i += 1
            line = lines[i]
        
        #print(in_geometry_block, in_bond_block, line)
        
        # parse geometry if appropriate
        if in_geometry_block:
            fields = re.split(' +', line)
            if len(fields)<6:
                print("Error parsing file:")
                print("Line = '%s'" % line.strip())
                print(fields)
                break
            x,y,z = float(fields[2]), float(fields[3]), float(fields[4])
            this_geometry.append([x,y,z])
            symbol = fields[5]
            this_symbol_list.append(symbol)
        elif in_bond_block:
            fields = re.split(' +', line.strip())
            if len(fields) == 4:
                # parse bonds, checking that the bonds are increasing
                try:
                    this_bond_number = int(fields[0])
                    atom1 = int(fields[1])
                    atom2 = int(fields[2])
                    bond_order = int(fields[3])
                    if this_bond_number != bond_number + 1:
                        raise ValueError("non-sequential bond number")
                    bond_number = this_bond_number
                    this_bond = (atom1, atom2, bond_order)
                    this_bonds.append(this_bond)
                except:
                    # assume we have left the bond block
                    in_geometry_block = False
                    in_bond_block = False
            else:
                # we have left the bond block
                in_geometry_block = False
                in_bond_block = False
        
        # go to next line
        i += 1
        
        # store geometry and reinitialize if appropriate
        if i == len(lines) or ( not in_geometry_block and not in_bond_block and len(this_geometry) > 0 ):
            all_geometries.append(np.array(this_geometry))
            all_symbol_lists.append(this_symbol_list)
            all_bonds.append(this_bonds)
            this_geometry = []
            this_symbol_list = []
            this_bonds = []
            
    # return result  
    all_geometries = np.array(all_geometries)
    result = Geometries(all_geometries, all_symbol_lists, all_bonds)
    print("%s read." % str(result))
    return result

### GJF file format reader
# reads a bunch of gjf files into a Geometries object
def read_conformers_from_gjfs(input_mask):
    input_filenames = list(glob(input_mask))
    input_filenames.sort()

    # initialize some lists
    all_geometries = []
    all_symbol_lists = []
    all_bonds = None
    this_geometry = []
    this_symbol_list = []

    for filename in input_filenames:
        this_geometry = []
        this_symbol_list = []
        with open(filename, 'r') as filehandle:
            lines = filehandle.read().splitlines()
        i = 0
        blanks = 0
        previous_line_was_blank = False
        found_geometry = False
        while i < len(lines):
            # read the current line
            line = lines[i].strip()

            # check for blank sections
            if len(line)==0 and not previous_line_was_blank:
                blanks += 1
                previous_line_was_blank = True
            else:
                previous_line_was_blank = False

            # skip to geometry section
            if blanks == 2 and not found_geometry:
                # skip charge and multiplicity card
                i += 2
                line = lines[i].strip()
                found_geometry = True
            if blanks == 3:
                break

            # parse geometry
            if blanks == 2:
                fields = re.split(' +', line)
                try:
                    x,y,z = float(fields[1]), float(fields[2]), float(fields[3])
                    this_geometry.append([x,y,z])
                    symbol = fields[0]
                    this_symbol_list.append(symbol)
                except:
                    print("Error parsing file:")
                    print("Line = '%s'" % line.strip())
                    print(fields)
                    raise ValueError()

            # go to next line
            i += 1

        all_geometries.append(np.array(this_geometry))
        all_symbol_lists.append(this_symbol_list)
    all_geometries = np.array(all_geometries)
    return Geometries(all_geometries,all_symbol_lists,all_bonds)

### Gaussian output file reader
def read_geometries_from_gaussian_out(input_mask, read_intermediate_geometries = False, successful_jobs_only = True):
    input_filenames = list(glob(input_mask))
    input_filenames.sort()

    # initialize some lists
    all_geometries = []
    all_symbol_lists = []
    all_bonds = None
    all_energies = []
    
    for filename in input_filenames:
        # initialize some arrays
        file_geometries = []
        file_symbol_lists = []
        file_energies = []

        this_geometry = []
        this_symbol_list = []
        this_energy = None

        with open(filename, 'r') as filehandle:
            lines = filehandle.read().splitlines()
        if successful_jobs_only and not lines[-1].strip().startswith("Normal termination"):
            print("Skipping %s as the job did not terminate normally." % filename)
            continue

        i = 0
        in_geometry_block = False
        while i < len(lines):
            # read the current line
            line = lines[i].strip()

            # detect geometry block
            if line == "Standard orientation:":
                i += 5
                in_geometry_block = True
                continue
            elif in_geometry_block and line.startswith("------"):
                i += 1
                in_geometry_block = False
                continue

            # read geometry if applicable
            if in_geometry_block:
                fields = re.split(' +', line)
                if len(fields) != 6:
                    print("error parsing >>> " + line)
                    raise ValueError("unexpected number of fields on geometry line")
                if len(this_geometry) > 0 and fields[0] == "1":
                    # reset fields
                    this_geometry = []
                    this_symbol_list = []
                    this_energy = None
                try:
                    x,y,z = float(fields[3]), float(fields[4]), float(fields[5])
                    this_geometry.append([x,y,z])
                    symbol = get_symbol(fields[1])
                    this_symbol_list.append(symbol)
                except:
                    print("error parsing >>> %s" % line)
                    print(fields)
                    raise ValueError("error parsing geometry")

            # read energy if applicable
            if not in_geometry_block and line.startswith("SCF Done"):
                fields = re.split(' +', line)
                if len(fields) != 9:
                    print("error parsing >>> " + line)
                    raise ValueError("unexpected number of fields on energy line")
                this_energy = float(fields[4])
            
                if len(this_geometry) == 0:
                    raise ValueError("energy without geometry")
                
                # store results
                #print(filename, len(this_geometry), len(this_symbol_list), this_energy, this_geometry[0])
                file_geometries.append(this_geometry)
                file_symbol_lists.append(this_symbol_list)
                file_energies.append(this_energy)

                # reinitialize arrays
                this_geometry = []
                this_symbol_list = []
                this_energy = None
 
            # go to next line
            i += 1

        # store results
        if read_intermediate_geometries:
            # store all geometries and energies from this file
            all_geometries.extend(file_geometries)
            all_symbol_lists.extend(file_symbol_lists)
            all_energies.extend(file_energies)
        elif len(file_geometries) > 0:
            # only store the last geometry and energy
            # (assuming there is something to store)
            all_geometries.append(file_geometries[-1])
            all_symbol_lists.append(file_symbol_lists[-1])
            all_energies.append(file_energies[-1])

        # reinitialize lsits
        file_geometries = []
        file_symbol_lists = []
        file_energies = []

    # return result
    return Geometries(all_geometries,all_symbol_lists,all_bonds,all_energies)

   


