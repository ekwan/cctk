"""
Functions to assist in optimizing structures.
"""

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
    if shutil.which(command) is not None:
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
        return run_xtb(molecule, nprocs=nprocs, return_energy=return_energy, opt=True)

def get_energy(molecule, method=Methods.GFN2_XTB, nprocs=1):
    """
    Dispatcher method to connect method to molecule.

    Args:
        molecule (cctk.Molecule):
        method (Methods):
        nprocs (int): number of cores to employ

    Returns:
        energy
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(method, Methods), "need a valid molecule!"

    if method is Methods.GFN2_XTB:
        _, energy = run_xtb(molecule, nprocs=nprocs, return_energy=True, opt=False)
        return energy

def run_xtb(molecule, nprocs=1, return_energy=False, opt=False):
    """
    Run ``xtb`` in a temporary directory and return the output molecule.
    """
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(nprocs, int)

    assert installed("xtb"), "xtb must be installed!"

    command = f"xtb --gfn 2 --chrg {molecule.charge} --uhf {molecule.multiplicity - 1}"
    if nprocs > 1:
        command += f" --parallel {nprocs}"

    if opt:
        command += " xtb-in.xyz --opt tight &> xtb-out.out"
    else:
        command += " xtb-in.xyz &> xtb-out.out"

    try:
        os.environ["OMP_NUM_THREADS"] = str(nprocs)
        os.environ["MKL_NUM_THREADS"] = str(nprocs)
        with tempfile.TemporaryDirectory() as tmpdir:
            cctk.XYZFile.write_molecule_to_file(f"{tmpdir}/xtb-in.xyz", molecule)
            sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=tmpdir, shell=True)

            output_mol, energy = None, None
            if opt:
                output_mol = cctk.XYZFile.read_file(f"{tmpdir}/xtbopt.xyz").get_molecule()
                energy_file = cctk.File.read_file(f"{tmpdir}/xtbopt.log")
                fields = energy_file[1].split()
                energy, gradient = float(fields[1]), float(fields[3])

            else:
                # stopgap solution but should work ok. XTB output files should actually be parsed eventually. 
                # ccw 4.15.21
                output_file = cctk.File.read_file(f"{tmpdir}/xtb-out.out")
                r = re.compile("total energy\s+(-?\d+.\d+)", re.IGNORECASE)
                for line in output_file[::-1]:
                    m = r.search(line)
                    if m:
                        energy = float(m.group(1))
                        break

            if return_energy:
                return output_mol, energy
            else:
                return output_mol
    except Exception as e:
        raise ValueError(f"Error running xtb:\n{e}")

def csearch(use_tempdir=True, **kwargs):
    """
    Run a conformational search on a molecule using ``crest``.

    Args:
        molecule (cctk.Molecule): molecule of interest
        constraints (list): list of atom numbers to constrain
        nprocs (int): number of processors to use
        noncovalent (Bool): whether or not to use non-covalent settings
        logfile (str): file to write ongoing ``crest`` output to
        additional_flags (str): flags to pass to command line

    Returns:
        cctk.ConformationalEnsemble
    """
    assert installed("crest"), "crest must be installed!"

    ensemble = None
    try:
        if use_tempdir:
            with tempfile.TemporaryDirectory() as tmpdir:
                ensemble = _do_csearch(directory=tmpdir, **kwargs)
        else:
            ensemble = _do_csearch(directory=os.getcwd(), **kwargs)
    except Exception as e:
        raise ValueError(f"Error running xtb:\n{e}")

    return ensemble

def _do_csearch(molecule, directory, gfn=2, nprocs=1, logfile=None, noncovalent=False, constraints=None, additional_flags=None):
    assert isinstance(molecule, cctk.Molecule), "need a valid molecule!"
    assert isinstance(nprocs, int)
    assert isinstance(logfile, str)

    assert gfn in [2, 1, "ff"], "invalid value for ``gfn``!"

    cctk.XYZFile.write_molecule_to_file(f"{directory}/xtb-in.xyz", molecule)

    nci = ""
    if noncovalent:
        nci = "-nci"

    command = None
    if constraints is not None:
        assert isinstance(constraints, list)
        assert all(isinstance(n, int) for n in constraints)
        command = f"crest xtb-in.xyz --constrain {','.join([str(c) for c in constraints])}"
        result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=directory, shell=True)
        result.check_returncode()
        command = f"crest xtb-in.xyz --gfn{gfn} --chrg {molecule.charge} -cinp .xcontrol.sample --uhf {molecule.multiplicity - 1} -T {nprocs} {nci}"
    else:
        command = f"crest xtb-in.xyz --gfn{gfn} --chrg {molecule.charge} --uhf {molecule.multiplicity - 1} -T {nprocs} {nci}"

    if additional_flags is not None:
        command = command + " " + additional_flags

    if logfile:
        with open(logfile, "w") as f:
            result = sp.run(command, stdout=f, stderr=f, cwd=directory, shell=True)
    else:
        result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=directory, shell=True)
    result.check_returncode()

    ensemble = cctk.XYZFile.read_ensemble(f"{directory}/crest_conformers.xyz").ensemble
    return ensemble



