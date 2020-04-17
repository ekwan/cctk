:orphan:

.. _about:

.. |br| raw:: html

  <br/>

.. |nbsp| raw:: html 

   &nbsp;

=====
About
=====

*cctk* simplifies routine tasks in computational quantum chemistry such as:

- file parsing and conversion (Gaussian ``.gjf`` and ``.out``, Orca ``.inp``, ``.mol2``, ``.xyz``, ``.pdb``, ``.mae``) 

- generating input files:
    - screening of methods and basis sets
    - automated resubmission of failed jobs
    - replacement of one functional group with another
    - systematic enumeration of conformers
    - potential energy surface scans

-  output file analysis:
    - job completion status and timings
    - extraction of geometries during optimizations or along IRCs
    - determination of bond lengths, angles, or dihedrals
    - calculating NMR shifts via linear scaling
    - molecular properties along reaction coordinates (e.g., charges)

*cctk* was written in Python 3 by Corin Wagen and Eugene Kwan to replace
an endless number of *ad hoc* scripts.

**How to cite:**

Wagen, C.C.; Kwan, E.E. *cctk* **2020**, `www.github.com/ekwan/cctk <https//www.github.com/ekwan/cctk>`_.

