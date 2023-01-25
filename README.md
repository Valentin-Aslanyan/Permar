# Permar - a suite of analysis routines for the Lare3d code

**LARE**S **PERMAR**INI: Roman guardian deities of seafarers

This package contains native Python3 routines to read, process and display outputs from the Lare3d code. **To make them work, you must change the path near the top of each file, which by default is `'/Change/This/Path'` to the directory where the file `Permar_Functions_Python.py` is kept.**

In particular, the SDF file format of the Centre for Fusion Space and Astrophysics (CFSA) at the University of Warwick (as opposed to a completely different and more widespread scientific SDF file format) can be read using a single Python function without the need for external routines.

There is also an extension to the Lare3d code itself to output the velocity at *e.g.* the bottom boundary every timestep, thereby allowing the “footpoints” of magnetic field lines to be accurately tracked throughout the simulation. Higher level quantities such as the squashing factor *Q* and the connectivity of magnetic field lines can also be calculated.

The package is written mostly in Python version ≥ 3.7, with additional high-performance routines in Fortran 90.

**Please acknowledge the use of these routines if you publish or present any results**. Feel free to get in touch, I am quite approachable.
