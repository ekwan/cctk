"""
Functions to assist in optimizing structures.
"""

import numpy as np
import os, tempfile, shutil, re
import cctk
import subprocess as sp

from enum import Enum

class Methods(Enum):
    """
    Enum of different computational methods. For now, just GFN2-xtb is implemented.
    """
    GFN2_XTB = "xtb"

def installed(command):
    if shutil.which("command") is not None:
        return True
    if re.search(command, os.environ["PATH"]):
        return True

    return False

def optimize_molecule(molecule, method=Methods.GFN2_XTB, nprocs=1, return_energy=False):
    """
    Dispatcher method to connect method to molecule.

    Args:
        molecule (cctk.Molecule):
        method (Methods):
        nprocs (int): number of cores to employ
        return_energy (Bool): to return energy or not

    Returns:
        molecule
        energy (optional)
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(method, Methods), "need a valid molecule!"

    if method is Methods.GFN2_XTB:
        return run_xtb(molecule, nprocs=nprocs, return_energy=return_energy)

def run_xtb(molecule, nprocs=1, return_energy=False):
    """
    Run ``xtb`` in a temporary directory and return the output molecule.
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(nprocs, int)

    assert installed("xtb"), "xtb must be installed!"

    command = f"xtb --gfn 2 --chrg {molecule.charge} --uhf {molecule.multiplicity - 1}"
    if nprocs > 1:
        command += f" --parallel {nprocs}"
    command += " xtb-in.xyz --opt tight &> xtb-out.out"

    output_mol, energy, gradient = None, None, None
    try:
        os.environ["OMP_NUM_THREADS"] = str(nprocs)
        os.environ["MKL_NUM_THREADS"] = str(nprocs)
        with tempfile.TemporaryDirectory() as tmpdir:
            cctk.XYZFile.write_molecule_to_file(f"{tmpdir}/xtb-in.xyz", molecule)
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=tmpdir, shell=True)

            output_mol = cctk.XYZFile.read_file(f"{tmpdir}/xtbopt.xyz").molecule
            energy_file = cctk.File.read_file(f"{tmpdir}/xtbopt.log")

            fields = energy_file[1].split()
            energy, gradient = float(fields[1]), float(fields[3])
    except Exception as e:
        raise ValueError(f"Error running xtb:\n{e}")

    if return_energy:
        return output_mol, energy
    else:
        return output_mol

def csearch(molecule, constrain_atoms=None, rotamers=True, nprocs=1):
    """
    Run a conformational search on a molecule using ``crest``.

    Args:
        molecule (cctk.Molecule): molecule of interest
        constraints (list): list of atom numbers to constrain
        rotamers (bool): return all rotamers or group into distinct conformers
        nprocs (int): number of processors to use

    Returns:
        cctk.ConformationalEnsemble
    """

    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(nprocs, int)

    assert installed("crest"), "crest must be installed!"

    ensemble = None
    try:
        if 1:
            tmpdir = "/n/home03/cwagen/sw/cctk/crest/"
            if constrain_atoms is not None:
                assert isinstance(constrain_atoms, list)
                assert all(isinstance(n, int) for n in constrain_atoms)
                command = f"crest --constrain {','.join(constrain_atoms)}"
                result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=tmpdir, shell=True)
                result.check_returncode()

            cctk.XYZFile.write_molecule_to_file(f"{tmpdir}/xtb-in.xyz", molecule)
            command = f"crest xtb-in.xyz --chrg {molecule.charge} --uhf {molecule.multiplicity - 1} -T {nprocs}"
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=tmpdir, shell=True)
            result.check_returncode()

            if rotamers:
                ensemble = cctk.XYZFile.read_ensemble(f"{tmpdir}/crest_rotamers.xyz")
            else:
                ensemble = cctk.XYZFile.read_ensemble(f"{tmpdir}/crest_conformers.xyz")

    except Exception as e:
        raise ValueError(f"Error running xtb:\n{e}")

    return ensemble

