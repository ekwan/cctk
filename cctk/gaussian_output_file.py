import sys
import re
import numpy as np

from cctk import OutputFile
from cctk.helper_functions import get_symbol, search_for_block, compute_distance_between, compute_angle_between, compute_dihedral_between, get_number

class GaussianOutputFile(OutputFile):
    '''
    Generic class for all Gaussian output files. 

    Attributes: 
        title (str): the title from the Gaussian file
        theory (dict): contains information from header
        header (str):  header from input file (e.g. "#p opt freq=noraman b3lyp/6-31g(d)").
        footer (str):  footer from input file, if any. genecp commands or opt=modredundant specifications go here. 
        success (Bool): if the file terminated successfully
    '''
    def __init__(self, filename):
        self.successful = False
        self.read_file(filename)    

    def read_file(self, filename):
        '''
        Read a Gaussian output file.
        Automatically determines the theory line and if the job was successful. 
        Returns an array of geometries and energies. 
        '''        
        lines = super().read_file(filename)

        file_geometries = []
        file_symbol_lists = []
        file_energies = []

        self.header = search_for_block(lines, "#p", "----")

        if lines[-1].strip().startswith("Normal termination"):
            self.successful = True

        return lines

    def read_geometries_and_energies(self, lines):
        """
        Reads geometries, symbol lists, and energies from the file.
        
        Args:
            lines (list): the list of lines in the file

        Returns:
            array of geometries (each of which itself is an array of arrays)
            list of atomic symbols
            array of energies
            array containing the number of SCF iterations per step 
        """
        
        file_geometries = []
        file_symbol_lists = []
        file_energies = []
        file_scf_iterations = []

        this_geometry = []
        this_symbol_list = []
        this_energy = None

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
                num_cycles = int(fields[7])
                 
                if len(this_geometry) == 0:
                    raise ValueError("energy without geometry")
                
                # store results
                file_geometries.append(this_geometry)
                file_symbol_lists.append(this_symbol_list)
                file_energies.append(this_energy)
                file_scf_iterations.append(num_cycles)

                # reinitialize arrays
                this_geometry = []
                this_symbol_list = []
                this_energy = None
 
            # go to next line
            i += 1

        # return result
        return file_geometries, file_symbol_lists[0], file_energies, file_scf_iterations

    def read_bonds(self, lines):
        """
        Reads bonding information from the output file. 
        
        Args:
            lines (list): the list of lines in the file

        Returns:
        """
        
        bond_array = []
        current_array = []

        i = 0 
        in_bonding_section = False
        while i < len(lines):
            # read the current line
            line = lines[i].strip()
            
            if in_bonding_section == False:
                if re.search(r"Initial Parameters", line):
                    i += 5
                    in_bonding_section = True
                    continue
                else: 
                    i += 1
                    continue
            else:
                if re.search(r"! R", line):
                    pieces = list(filter(None, line.split(" ")))
                    atoms = pieces[2].replace("R", '').replace("(", "").replace(")","").split(",")
                    try:
                        current_array.append([int(atoms[0]), int(atoms[1])])
                    except:
                        raise ValueError(f"error parsing line {i} - can't extract atoms!")
                    i += 1
                    continue
                elif re.search(r"! A", line):
                    in_bonding_section = False
                    bond_array = current_array
                    current_array = []
                    i += 1
                    continue
                else: 
                    raise ValueError(f"can't parse line {i} for bonding!")

        return bond_array

    def _find_parameter(self, lines, parameter, expected_length, which_field):
        """
        Helper method to search through the output file and find key forces and displacements. 

        Args:
            lines (list): list of lines in file
            parameter (string): test to search for
            expected_length (int): how many fields there should be
            which_field (int): which field the parameter is (zero-indexed) 
        Returns:
            a list of all the extracted values
        """
        
        if (not isinstance(expected_length, int)) or (not isinstance(which_field, int)):
            raise TypeError("expected_length and which_field must be type int!")
     
        if which_field >= expected_length:
            raise ValueError("can't expect a field after the last field!") 
       
        matches = []
        pattern = False
        
        try:
            pattern = re.compile(parameter)
        except: 
            raise ValueError("pattern {pattern} cannot be compiled as a regex; try again!")
        
        if pattern:
            for line in lines:
                if pattern.search(line):
                    fields = re.split(' +', line)
                    fields = list(filter(None, fields))
                 
                    if len(fields) == expected_length:
                        matches.append(float(fields[which_field]))
            return matches

    def get_distance(self, geometry, atom1, atom2):
        """
        Wrapper to compute distance between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
        except:
            raise TypeError("atom1 and atom2 must be int!")
        
        if (atom1 <= len(geometry)) and (atom2 <= len(geometry)):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            return compute_distance_between(v1, v2)
        else:
            raise ValueError(f"atom numbers {atom1} or {atom2} too big!")

    def get_angle(self, geometry, atom1, atom2, atom3):
        """
        Wrapper to compute angle between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
            atom3 = int(atom3)
        except:
            raise TypeError("atom1, atom2, and atom3 must be int!")
        
        if ((atom1 <= len(geometry)) and (atom2 <= len(geometry))) and (atom3 <= len(geometry)):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            v3 = np.array(geometry[atom3-1])
            return compute_angle_between(v1, v2, v3)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, or {atom3} too big!")
   
    def get_dihedral(self, geometry, atom1, atom2, atom3, atom4):
        """
        Wrapper to compute dihedral angle between two atoms. 
        """
        try:
            atom1 = int(atom1)
            atom2 = int(atom2)
            atom3 = int(atom3)
            atom4 = int(atom4)
        except:
            raise TypeError("atom1, atom2, atom3, and atom4 must be int!")
        
        if ((atom1 <= len(geometry)) and (atom2 <= len(geometry))) and ((atom3 <= len(geometry)) and (atom4 <= len(geometry))):
            v1 = np.array(geometry[atom1-1])
            v2 = np.array(geometry[atom2-1])
            v3 = np.array(geometry[atom3-1])
            v4 = np.array(geometry[atom4-1])
            return compute_dihedral_between(v1, v2, v3, v4)
        else:
            raise ValueError(f"atom numbers {atom1}, {atom2}, {atom3}, or {atom4} too big!")
