.. _tutorial_01:

================================
Tutorial 01: Changing File Types
================================

Objectives
==========

This tutorial will teach:

- Creating command-line *cctk* scripts.
- Manipulation of ``File`` objects.
- Reading/writing data.

Overview
========

Many computational chemistry papers report structures in the ``.xyz`` format, which is not recognized by Gaussian. 
Although manual conversion is facile for one file, an automated solution can prove useful for bulk data processing.
A sample ``.xyz`` file can be found at the end of this tutorial. 

This tutorial will showcase *cctk*'s ability to automatically interconvert between file types, as well as provide a template for other command-line scripts.

Creating a Bash Script
======================

Open a terminal window in a directory with a file titled ``tutorial1.xyz``. A sample xyz file can be found here `tutorial1.xyz <../tutorial/tutorial_01/tutorial1.xyz>`_.
In a terminal window, create a new file called ``read_from_xyz_01.py`` and open it in your favorite text editor (e.g., ``vim``, ``emacs``, or ``nano``).
::

    $ vim read_from_xyz_01.py

This will open a blank file. First, we need to load *cctk*::

    import re
    from cctk import XYZFile, GaussianFile

Now that we've loaded *cctk*, we can read in data from an input file::

    filename = "./tutorial1.xyz"
    file = XYZFile.read_file(filename)

The above code creates a *cctk* ``XYZFile`` object from the file we specified, which now exists as a Python data structure. 
To output a different filetype, we need to extract the ``Molecule`` object represented by the file and write it as a ``.gjf`` file::

    molecule = file.get_molecule()

    newfile = filename.rsplit('/',1)[-1]
    newfile = re.sub(r"xyz$", "gjf", newfile)

    GaussianFile.write_molecule_to_file(
        newfile,
        molecule,
        "#p opt freq=noraman b3lyp/6-31g(d) empiricaldispersion=gd3bj",
        None,
    )

The command ``write_molecule_to_file`` is a class method, meaning we can create a ``.gjf`` file without needing to create another Python object. 
All we need to supply is the path to the new file, the ``Molecule`` object (in this case, ``file.molecule``), and the header and footer for the new file. 
(In this case, we have also ensured that the output file ends in ``.gjf`` and is placed in the directory from which we run the script by using Python string manipulation.)

Running this script on `tutorial1.xyz <./../test/static/tutorial1.xyz>`_ generates the desired ``test.gjf`` input file::

    $ python read_from_xyz.py

The start of `tutorial1.gjf <./tutorial1.gjf>`_ is shown below::

    %nprocshared=16
    %mem=32GB
    #p opt freq=noraman b3lyp/6-31g(d) empiricaldispersion=gd3bj

    title

    0 1
    6 0.25892000 0.68427000 0.00004500
    ...

The and the start of `tutorial1.inp <./tutorial1.inp>`_ ::

    ! opt freq b3lyp/6-31g(d) d3bj
    %maxcore 2000
    %pal
	    nproc 16
    end

    * xyz 0 1
    6          0.25892001    0.68427002    0.00004500
    ...

The script works!

Adding Command-Line Arguments
=============================

To create a more user-friendly script, we might want to make it so that we can specify the file and desired header without manually editing the script each time. 
This can be done using Python's ``argparse`` module::

    import sys, argparse, re
    from cctk import GaussianFile, XYZFile
    
    parser = argparse.ArgumentParser(prog="resubmit.py")
    parser.add_argument("--header", "-h", type=str)
    parser.add_argument("filename")
    args = vars(parser.parse_args(sys.argv[1:]))

    assert args["filename"], "Can't read file without a filename!"
    assert args["header"], "Can't write file without a header!"

The script will now expect two arguments, the first of which must be preceded by the ``-h`` flag. 

After adding comments and integrating the above variables throughout, the final script looks like this::

    import sys, argparse, re
    from cctk import GaussianFile, XYZFile

    #### Usage: python read_from_xyz.py -h "#p opt freq=noraman b3lyp/6-31g(d)" path/to/file.xyz

    parser = argparse.ArgumentParser(prog="resubmit.py")
    parser.add_argument("--header", "-h", type=str)
    parser.add_argument("filename")
    args = vars(parser.parse_args(sys.argv[1:]))

    assert args["filename"], "Can't read file without a filename!"
    assert args["header"], "Can't write file without a header!"

    file = XYZFile.read_file(args["filename"])
    newfile = args["filename"].rsplit('/',1)[-1]
    newfile = re.sub(r"xyz$", "gjf", newfile)

    GaussianFile.write_molecule_to_file(
        newfile,
        file.molecule,
        args["header"],
        None,
    )

To run this on our test file, simply type::

    python read_from_xyz.py -h "#p opt b3lyp/6-31(g)" test.xyz

This script can now be copied to other directories and used as a command-line tool.
The template provided here can also be modified for myriad *cctk*-based applications, as future tutorials will demonstrate.
