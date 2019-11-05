import sys
import re
import numpy as np

from cctk import InputFile

class GaussianInputFile(InputFile):
    """
    Class for Gaussian input files 

    Attributes:
        atoms (list): list of atomic numbers 
        geometry (list): list of geometries for each cycle (so really a list of lists)
        theory (dict): optional, dictionary of values related to level of theory
        header (str): optional, header of .gjf file
        footer (str): optional, footer of .gjf file
        title (str): optional, title of .gjf file
        charge (int): the charge of the molecule
        multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)
    """

    def __init__(self, atoms, geometry, theory=None, header=None, footer=None, title='title', charge=0, multiplicity=1):
        """
        Create new GaussianInputFile object.
        """
        if len(atoms) != len(geometry):
            raise ValueError("length of atoms and length of geometry arrays must be the same!")

        self.atoms = [int(z) for z in atoms]
        self.geometry = geometry
        self.theory = theory
        self.header = header
        self.footer = footer
        self.title = title
        self.charge = charge
        self.multiplicity = multiplicity

    def write_file(self, filename, memory=32, cores=16, chk_path=None):
        """ 
        Write a .gjf file, using object attributes. 
        """
        
        if self.header:
            text = ''
            text += "%nprocshared={}GB\n".format(cores) 
            text += "%mem={}GB\n".format(memory) 
            
            if chk_path:
                text += "%chk={}\n".format(chk_path)   
            
            text += self.header.rstrip()
            text += "\n"
            text += "\n"
            text += "{}\n".format(self.title)
            text += "\n"
           
            text += "{} {}\n".format(self.charge, self.multiplicity)
            for index, line in enumerate(self.geometry):
                text += "{:2d} {:.8f} {:.8f} {:.8f}\n".format(self.atoms[index], line[0], line[1], line[2])
            
            text += "\n" 
            if self.footer:
                text += self.footer.rstrip()
                text += "\n"
                text += "\n"

            super().write_file(filename, text)
        else: 
            raise ValueError("need header to write input file!")
