import numpy as np

from cctk import GaussianInputFile

class GaussianJob():
    """
    Creates input files of the specific type through factory methods.
    """
    
    @staticmethod
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

    @staticmethod
    def create_opt(cls, atoms, geometry, functional='b3lyp', basis='6-31g(d)', freq=True, scrf=False, solvent='', disp=False, charge=0, multiplicity=1):
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
            disp (str): if False, no empirical dispersion will be used. if a string, then that string will be used as the empirical dispersion keyword (e.g. ``gd3bj``). 
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        Returns:
            the created GaussianInputFile object
        """
        cls._check_atoms_geometry(atoms, geometry)
        
        header = "#p opt"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += f" scrf=({scrf}, solvent={solvent})" 
        if disp:
            header += f" empiricaldispersion={disp}"
        
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=None, charge=charge, multiplicity=multiplicity)
    
    @staticmethod
    def create_ts_opt(cls, atoms, geometry, functional='b3lyp', basis='6-31g(d)', freq=True, scrf=False, solvent='', disp=False, charge=0, multiplicity=1):
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
            disp (str): if False, no empirical dispersion will be used. if a string, then that string will be used as the empirical dispersion keyword (e.g. ``gd3bj``). 
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        Returns 
            the created GaussianInputFile object
        """
        cls._check_atoms_geometry(atoms, geometry)
        
        header = "#p opt=(ts, calcfc, noeigentest)"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += f" scrf=({scrf}, solvent={solvent})" 
        if disp:
            header += f" empiricaldispersion={disp}"

        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=None, charge=charge, multiplicity=multiplicity)

    @staticmethod
    def create_constrained_opt(cls, atoms, geometry, constraints, functional='b3lyp', basis='6-31g(d)', freq=False, scrf=False, solvent='', disp=False, charge=0, multiplicity=1):
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
            disp (str): if False, no empirical dispersion will be used. if a string, then that string will be used as the empirical dispersion keyword (e.g. ``gd3bj``). 
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        Returns 
            the created GaussianInputFile object
        """
        cls._check_atoms_geometry(atoms, geometry)
        if len(constraints) == 0:
            raise TypeError("no constraints -- perhaps a regular opt would be better?")
        
        header = "#p opt=modredundant"
        if freq:
            header += " freq=noraman"
        header += f" {functional}/{basis}"
        if (scrf and solvent):
            header += f" scrf=({scrf}, solvent={solvent})" 
        if disp:
            header += f" empiricaldispersion={disp}"

        footer = ''
        for constraint in constraints:
            footer += " ".join(str(x) for x in constraint)
            footer += "\n"    
        
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=footer, charge=charge, multiplicity=multiplicity)
    
    @staticmethod
    def create_nmr(cls, atoms, geometry, nmr="giao" functional='b3lyp', basis='6-31g(d)', scrf=False, solvent='', disp=False, charge=0, multiplicity=1):
        """
        Creates a Gaussian NMR job. 

        Args:
            atoms (list): list of atomic numbers 
            geometry (list): list of geometries for each cycle (so really a list of lists)
            nmr (str): the method by which the shielding tensors will be computed (e..g ``giao`` or ``csgt``)
            functional (str): the functional to use for the optimization
            basis (str): the basis set to use for the optimization
            scrf (str): if False, no implicit solvent calculation will be performed. if a string, that string will be used as the type of implicit solvation (e.g. ``smd`` or ``pcm``). 
            solvent (str): the solvent to use for implicit solvent calculation (refer to Gaussian notes for full list).  
            disp (str): if False, no empirical dispersion will be used. if a string, then that string will be used as the empirical dispersion keyword (e.g. ``gd3bj``). 
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        Returns:
            the created GaussianInputFile object
        """
        cls._check_atoms_geometry(atoms, geometry)
        if (not nmr) or (not isinstance(nmr, str)):
            raise TypeError("can't run nmr job without keyword nmr!")
        
        header = f"#p nmr={nmr} {functional}/{basis}"
        if (scrf and solvent):
            header += f" scrf=({scrf}, solvent={solvent})" 
        if disp:
            header += f" empiricaldispersion={disp}"
        
        return GaussianInputFile(atoms, geometry, theory=None, header=header, footer=None, charge=charge, multiplicity=multiplicity)
    
    @staticmethod
    def create_irc(cls, atoms, geometry, functional='b3lyp', basis='6-31g(d)', scrf=False, solvent='', disp=False, charge=0, multiplicity=1):
        """
        Creates two Gaussian IRC jobs, one in the forward direction and one in the reverse direction. The ``CalcFc`` option is used. 

        Args:
            atoms (list): list of atomic numbers 
            geometry (list): list of geometries for each cycle (so really a list of lists)
            functional (str): the functional to use for the optimization
            basis (str): the basis set to use for the optimization
            scrf (str): if False, no implicit solvent calculation will be performed. if a string, that string will be used as the type of implicit solvation (e.g. ``smd`` or ``pcm``). 
            solvent (str): the solvent to use for implicit solvent calculation (refer to Gaussian notes for full list).  
            disp (str): if False, no empirical dispersion will be used. if a string, then that string will be used as the empirical dispersion keyword (e.g. ``gd3bj``). 
            charge (int): the charge of the molecule
            multiplicity (int): the spin state of the molecule (1 corresponds to singlet, 2 to doublet, 3 to triplet, etc. -- so a multiplicity of 1 is equivalent to S=0)

        Returns:
            the created GaussianInputFile objects, forward first and then reverse
        """
        cls._check_atoms_geometry(atoms, geometry)
        
        header_f = f"#p irc=(calcfc, forward) {functional}/{basis}"
        header_r = f"#p irc=(calcfc, reverse) {functional}/{basis}"
        if (scrf and solvent):
            header_f += f" scrf=({scrf}, solvent={solvent})" 
            header_r += f" scrf=({scrf}, solvent={solvent})" 
        if disp:
            header_f += f" empiricaldispersion={disp}"
            header_r += f" scrf=({scrf}, solvent={solvent})" 
        
        job_f = GaussianInputFile(atoms, geometry, theory=None, header=header_f, footer=None, charge=charge, multiplicity=multiplicity)
        job_r = GaussianInputFile(atoms, geometry, theory=None, header=header_r, footer=None, charge=charge, multiplicity=multiplicity)
        return job_f, job_r
