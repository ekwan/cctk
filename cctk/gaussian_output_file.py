import sys
import re
import numpy as np

from cctk import OutputFile
from cctk.helper_functions import get_symbol, search_for_block

class GaussianOutputFile(OutputFile):
    '''
    title = title from file
    theory = dict, containing information from header
    header = header from input file
    footer = footer from input file
    success = Bool, if 
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
        '''
        Reads geometries, symbol lists, and energies from the file.
        Returns an array of geometries, a SINGLE array of the symbols, an array of energies, and an array of # SCF cycles from the file.
        '''
        
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

    def geometry():
        pass
        
    def energy():
        pass

