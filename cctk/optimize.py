"""
Functions to assist in optimizing structures.
"""

import numpy as np
import os, tempfile, shutil
import cctk
import subprocess as sp

from enum import Enum

class Methods(Enum):
    """
    Enum of different computational methods. For now, just GFN2-xtb is implemented.
    """
    GFN2_XTB = "xtb"

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

    assert shutil.which("xtb") is not None, "xtb must be installed for optimization to work!"

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

