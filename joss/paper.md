
---
title: 'qtools: A Python toolset for the Q molecular simulation package.'  
tags:
  - Python
  - computational enzymology
  - empirical valence bond
  - EVB
  - Q
  - free energy perturbation
  - FEP  
authors:
  - name: Miha Purg
    orcid: 0000-0003-4647-6103
    affiliation: 1  
affiliations:
 - name: "Department of Chemistry-BMC, Uppsala University, Box 576, SE-751 23 Uppsala, Sweden"
   index: 1  
date: 4 September 2018  
bibliography: paper.bib
---
 
# Introduction

Molecular modeling and simulations are becoming an increasingly integral part of enzymology, complementing experimental measurements with a direct view into the dynamics and reactivity of enzymes at the atomic level. The multiscale quantum mechanics/molecular mechanics (QM/MM) methodology has proven to be particularly valuable for studying enzyme-catalyzed reactions, gaining widespread recognition evidenced by the Nobel prize in Chemistry in 2013.
By providing a map of the energy landscape for (bio)chemical systems, these methods help to elucidate chemical mechanisms, rationalize substrate specificity, and estimate energetic contributions to catalysis, all of which might otherwise be difficult to quantify.
Such detailed understanding of the forces and interactions governing the activity of enzymes is also tremendously valuable in the rational design of inhibitors as potential drug candidates. Furthermore, it can aid in the design of enzyme functionality, *e.g.* in search of environmentally friendly industrial biocatalysts.

One of the well-established methodologies for studying enzyme reactions is the empirical valence bond (EVB) method. [@warshel:1980] 
In EVB, a chemical reaction is described within a valence bond formalism, representing individual valence states (*e.g.* reactants and products) with classical force-field potentials. 
Employing a classical model dramatically reduces the computational cost compared to traditional QM/MM methods allowing for substantially longer molecular dynamics simulations, although at the cost of accuracy. A particularly useful application of EVB is in the rapid screening of point-mutations, as demonstrated in our recent work. [@amrein:2017] 
The EVB methodology is implemented in Q, a software package designed specifically for running molecular dynamics simulations in the context of free energy perturbation (FEP) and EVB. [@marelius:1998],[@bauer:2018] Due to its sleek and minimalistic design, the users would benefit greatly from a well-designed high-level interface, increasing the productivity and minimizing the use of system-specific, poorly tested, closed-source scripts.


# Overview

`qtools` is a toolset comprised of a Python library `Qpyl` and a set of command-line tools, for use with the molecular simulation package Q. The motivation behind `qtools` is to provide users of Q with a structured, well-tested, open-source Python library for efficient setup, troubleshooting, calibration, and analysis of EVB simulations. In addition, a set of Unix command-line tools, which directly calls the `Qpyl` API, is bundled in the repository, providing a user interface which is particularly practical for use in typical high-performance-computing (HPC) situations.
While EVB simulations are the main focus of `qtools`, non-reactive FEP calculations (*e.g.* ligand binding and solvation free energies) are fully supported.

The functionality implemented in `Qpyl` is application-specific but relatively broad, encompassing several modules for wrapping Q executables, generating and parsing I/O files, generating, parsing and converting parameter files, analysing FEP, EVB, QCP free energy results, calibrating the EVB potential, and other. 
Most importantly, it exposes this functionality *via* a Python API, providing users with programmable means for preparing, simulating, and analysing EVB simulations.
In addition to the EVB-related functionality, `qtools` implements the linear response approximation (LRA) method for estimating free energies. [@lee:1992] The LRA method defines the free energy between two states, A and B, as: 

$$
\Delta G(A \to B) = \frac{1}{2} \Big[ \langle U_B - U_A \rangle_A + \ \langle U_B - U_A \rangle_B \Big] 
$$

In `Qpyl`, LRA is used as an efficient and convenient post-processing approach for calculating the *virtual deletion* of amino-acid residues, providing estimates of their energetic contributions to catalysis. Additionally, `qpyl` allows the estimation of contributions to *reorganization energy*, [@muegge:1997] which were suggested to be especially useful in computational design of enzymes. [@fuxreiter:2014] 

Scientific publications employing `qtools` include [@amrein:2017],[@purg:2017],[@bauer:2018]. A recently published chapter in Methods in Enzymology provides step-by-step instructions on using `qtools` in practice, from obtaining force-field parameters and input generation, to analysis and troubleshooting. [@purg:2018]



# Technical details

The main design characteristics of `qtools` are:  

1. *Minimum dependencies*. The only dependency is Python 2.7, with the exception of one (non-essential) tool which requires `matplotlib`. 
2. *Logic separation and extensibility.* The functionality is implemented in the Python library `Qpyl` and exposed *via* an API. A higher-level application/interface (GUI, web-service, *etc.*) should avoid implementing Q-related functionality and instead call the API, analogous to the bundled command-line tools.
3. *Thorough testing*. Continuous integration using `py.test` for Qpyl, and Bash regression tests for command-line tools. 
4. *Strict sanity checks*. Invalid Q inputs, parameters and/or settings, charge imbalances, duplicate parameters, *etc.*, raise exceptions or propagate logging warnings. 
5. *Balanced abstraction*. `qtools` tries to maintain a good balance between the level of user-control over Q (*i.e.* relatively low abstraction level) and ease-of-use.

#### Python library Qpyl

Core functionality implemented in `Qpyl.core` includes generator and parser objects for Qfep, Qdyn, and Qcalc input and output files, as well as for Q library, parameter, and coordinate files.
Several external library and parameter files are supported, including AMBER [@maier:2015] `.parm`, `.frcmod`, `.lib` and `.prepin` files, and output from MacroModel. [@macromodel:2016]
`Qpyl.core.qtopology` implements an internal topology builder for mapping the system's bonding patterns with force-field parameters. In conjunction with the force-field potential energy functions implemented in `Qpyl.core.qpotential` one can use the topology object for direct evaluation of individual topological components of the system.  

###### Example 1: Joining parameter files
```python
>>> from Qpyl.core.qparameter import QPrm
>>> qprm = QPrm("oplsaa")
>>> qprm.read_prm("qoplsaa.prm")
>>> qprm.read_prm("ligand.prm") # raises QPrmError on parameter conflicts
>>> open("ooplsaa_ligand.prm", "w").write(qprm.get_string())
```

Higher-level functionality includes wrapper objects for Q programs, allowing for a thread-enabled automatic generation of inputs, system calls to Q executables, and analysis of outputs. 
The analysis objects in `Qpyl.qanalysis` extract/calculate relevant free energies (activation, reaction and FEP free energies) and perform basic statistical analysis of the simulations. Quantum classical path (QCP) and group exclusion calculations implemented in the recent version of Q are fully supported. 
In addition, `Qpyl` implements the automated calibration of (two-state) EVB potentials, the calculation of *group contributions* using the LRA method, and provides convenient methods for generating Qdyn input files for equilibration and FEP simulation stages.
`Qpyl.qmakefep` implements a semi-automated approach for determining the changes in force-field parameters between (valence) states, as defined in a Q FEP file.
The Q FEP file is built automatically from the provided parameter files (for each state), the initial state coordinate file, and a user-defined mapping of atoms in each state.
It further provides sanity checks for missing atoms and parameters, charge imbalances, *etc.*

###### Example 2: Analysing FEP outputs 
```python
>>> from Qpyl.qanalysis import QAnalyseFeps
>>> qaf = QAnalyseFeps(["1/qfep.out", "2/qfep.out"])
>>> print qaf.dgas.values()  # activation free energies
[12.21, 12.55]
```

###### Example 3: Generating a FEP file
```python
>>> from Qpyl.qmakefep import make_fep
>>> fepstring = make_fep(qmap_file="atom_mapping.qmap", pdb_file="sn2_start.pdb", forcefield="oplsaa",
                         parm_files=["cl.prm", "mebr.prm", "mecl.prm", "br.prm"],
                         lib_files=["cl.lib", "mebr.lib", "mecl.prm", "br.lib"])
>>> open("sn2.fep", "w").write(fepstring)
```

#### Command-line interface

`qtools` comes with a command-line interface allowing Q-users direct access to the functionality implemented in `Qpyl`. The CLI tools are designed for use in a Unix environment, are namespaced by using a prefix 'q_' (*e.g.* `q_analysefeps.py`), are non-interactive and rely extensively on command-line arguments.
Usages include: parameter conversion from MacroModel (OPLS) or AMBER formats (`q_ffld2q.py`, `q_amber2q.py`), FEP file generation (`q_makefep.py`), input generation (`q_genrelax.py`, `q_genfeps.py`), Qfep wrappers (`q_mapper.py`, `q_automapper.py`), Qcalc wrapper (`q_calc.py`), and analysis of simulations (`q_analysefeps.py`, `q_analysedyns.py`).

#### Graphical interface

Finally, `qtools` comes with a GUI application `q_plot.py` for visualizing results and troubleshooting simulations. `q_plot.py` displays the data generated by CLI tools, *i.e.* free energy profiles, LRA contributions, LRA group contributions, potential energy breakdown, kinetic energy, temperature, *etc.*, which is especially useful in troubleshooting simulations in Q. Some examples of use are demonstrated in reference [@purg:2018]


# Acknowledgments

I thank Dr. Janez Mavri and Dr. Paul Bauer for all the helpful discussions and support.

