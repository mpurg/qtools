"""This subpackage contains the core functionality of Qpyl.

Modules *qdyn*, *qfep*, and *qcalc* are used for wrapping Q executables,
generating inputs and parsing outputs.
Modules  *qparameter* and *qlibrary* are used for generating, converting,
and parsing Q parameter files. Reading of Amber and Macromodel formats is
also supported.
*qstructure* contains classes for reading and writing structure files.
*qtopology* implements an internal topology builder for mapping the system's bonding
patterns with the parameters. In conjunction with the force-field potential
functions implemented in *qpotential* one can use the topology object for direct
evaluation of individual topological components of the system.
"""

