import numpy as np

from cctk import GaussianInputFile

class GaussianJob():
    """
    Creates input files of the specific type through factory methods.
    """
    
    @classmethod
    def _check_atoms_geometry(cls, atoms, geometry):
        """
        Performs basic checks on ``atoms``and ``geometry`` prior to creating a ``GaussianInputFile``.   
        """
        if len(atoms) != len(geometry):
            raise ValueError("length of atoms and length of geometry arrays must be the same!")
       
        for atom in atoms: 
            if not isinstance(atom, int):
                raise TypeError(f"atom {atom} is not int")    

            if atom <= 0:
                raise ValueError(f"{atom} is not a valid atom number")

        for vector in geometry: 
            if len(vector) != 3:
                raise TypeError("each element of geometry must be a 3-element list!")
    @classmethod
    def create_opt(cls, atoms, geometry, functional='b3lyp', basis='6-31g(d)', freq=True, scrf=False, solvent='', charge=0, multiplicity=1):
        """
        Creates a basic Gaussian optimization job. 

        Args:
            atoms (list): list of atomic numbers 
            geometry (list): list of geometries for each cycle (so really a list of lists)
            functional (str): the functional to use for the optimization
            basis (str): the basis set to use for the optimization
            freq (Bool): whether to perform a subsequent frequency calculation or not
            scrf (str): if False, no implicit solvent calculation will be performed. if a string, that string will be used as the type of implicit solvation (e.g. ``smd`` or ``pcm``). 
            solvent (str): the solvent to use for implicit solvent calculation (refer to Gaussian notes for full list).  
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        """
        cls._check_atoms_geometry(atoms, geometry)
        
        header = "#p opt"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += "scrf=({scrf}, solvent={solvent})" 
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=None, charge=charge, multiplicity=multiplicity)
    
    @classmethod
    def create_ts_opt(cls, atoms, geometry, functional='b3lyp', basis='6-31g(d)', freq=True, scrf=False, solvent='', charge=0, multiplicity=1):
        """
        Creates a basic Gaussian transition state optimization job. 

        Args:
            atoms (list): list of atomic numbers 
            geometry (list): list of geometries for each cycle (so really a list of lists)
            functional (str): the functional to use for the optimization
            basis (str): the basis set to use for the optimization
            freq (Bool): whether to perform a subsequent frequency calculation or not
            scrf (str): if False, no implicit solvent calculation will be performed. if a string, that string will be used as the type of implicit solvation (e.g. ``smd`` or ``pcm``). 
            solvent (str): the solvent to use for implicit solvent calculation (refer to Gaussian notes for full list).  
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        """
        cls._check_atoms_geometry(atoms, geometry)
        
        header = "#p opt=(ts, calcfc, noeigentest)"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += "scrf=({scrf}, solvent={solvent})" 
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=None, charge=charge, multiplicity=multiplicity)

    @classmethod
    def create_constrained_opt(cls, atoms, geometry, constraints, functional='b3lyp', basis='6-31g(d)', freq=False, scrf=False, solvent='', charge=0, multiplicity=1):
        """
        Creates a basic Gaussian optimization job. 

        Args:
            atoms (list): list of atomic numbers 
            geometry (list): list of geometries for each cycle (so really a list of lists)
            constraints (list): list of constraints to employ (e.g. ``[['B', 5, 10, 'F'], ['B', 11, 12, 'S', 100, 0.01]]``)
            functional (str): the functional to use for the optimization
            basis (str): the basis set to use for the optimization
            freq (Bool): whether to perform a subsequent frequency calculation or not
            scrf (str): if False, no implicit solvent calculation will be performed. if a string, that string will be used as the type of implicit solvation (e.g. ``smd`` or ``pcm``). 
            solvent (str): the solvent to use for implicit solvent calculation (refer to Gaussian notes for full list).  
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        """
        cls._check_atoms_geometry(atoms, geometry)
        if len(constraints) == 0:
            raise TypeError("no constraints -- perhaps a regular opt would be better?")
        
        header = "#p opt=modredundant"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += "scrf=({scrf}, solvent={solvent})" 

        footer = ''
        for constraint in constraints:
            footer += " ".join(str(x) for x in constraint)
            footer += "\n"    
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=footer, charge=charge, multiplicity=multiplicity)
