### qtools 

[![Travis](https://api.travis-ci.org/mpurg/qtools.svg?branch=master)](https://travis-ci.org/mpurg/qtools)  [![Coverage Status](https://coveralls.io/repos/github/mpurg/qtools/badge.svg?branch=master)](https://coveralls.io/github/mpurg/qtools?branch=master)

#### Description

qtools is a toolset comprised of a Python library Qpyl and a set of command-line tools, for use with the [molecular dynamics simulation package *Q*](http://xray.bmc.uu.se/~aqwww/q/). This open-source project was created by Miha Purg (currently at Uppsala University) in order to provide Q-users with a set of well-tested tools for automating repetitive and/or tedious tasks, debugging simulations, and analysing results. The author strives to achieve a good balance between the level of user-control over *Q* (i.e., relatively low abstraction level) and the ease-of-use provided by qtools.

*Disclaimer: qtools is not affiliated with *Q* nor should it be considered an official toolset for Q.*


#### Functionality

- Python wrapper for Q
- Parameter conversions (FFLD OPLS, Amber/GAFF)
- Input generation for molecular dynamics/free energy perturbation/umbrella sampling (MD/FEP/US) simulations
- FEP-file generation for empirical valence bond (EVB) simulations
- EVB mapping and automatic fitting of off-diagonal and gas-shift parameters
- Linear response approximation (LRA) calculations (including group contributions)
- Output parsing and analysis of MD, FEP, and EVB
- Support for Q6 (quantum classical path and group exclusion calculations)
- Graphical tool for visualizing results and debugging simulations:
  - Free-energy profiles & sampling distributions
  - LRA and reorganization free-energy breakdown
  - LRA group contributions
  - Energy breakdown (vdW, electrostatics, ...)


#### Requirements

The library and tools do not have any dependencies, except *matplotlib*
which is required only for *q_plot.py*, a non-essential utility used mostly for
troubleshooting simulations.

The only requirements are a *Unix-like OS (Linux or OSX)*, *Python (versions 2.7, 3.4+)* and *Q* (see the [Q6 github page](https://github.com/qusers/Q6)).

#### Installation (Linux, OSX)

Clone the repository to your local directory:  
```
mkdir -p ~/bin && cd ~/bin
git clone https://github.com/mpurg/qtools
```

Add this line to your `~/.bashrc` or `~/.bash_profile`:  
```
source $HOME/bin/qtools/qtools_init.sh
```

Open a new shell session and run the CLI config script:  
```
qscripts_config.py
```
This script looks for Q executables in your `PATH` and creates a default settings-file `$QTOOLS_HOME/qscripts.cfg`.  
If the executables are not found the settings file should be modified manually.

#### Documentation

The starting point should be the `docs/tutorials` folder.
It contains working examples of how you can use the command line tools to prepare, simulate and analyse a simple SN2 reaction with EVB. Additionally, a recently published chapter in the serial Methods in Enzymology provides step-by-step instructions for performing EVB simulations on a real enzymatic system with `qtools`. Manuscript available [here](https://doi.org/10.1016/bs.mie.2018.06.007) (paywalled).
  
For those inclined on using the Python API directly, jupyter-notebook examples calling the Qpyl API can be found in `docs/examples`. The API reference can be found in `docs/api/html`. Additionally, the tests, and the command-line interface (qscripts-cli) can be useful to explore the API.
  
#### Testing

Automated Testing is performed on [Travis-CI](https://travis-ci.org/mpurg/qtools),
with [pytest](https://docs.pytest.org/en/latest/) for Qpyl
([code coverage](https://coveralls.io/github/mpurg/qtools?branch=master)).
and simple regression tests for CLI tools.
To run the tests locally make sure you have `pytest` and `pytest-cov` installed, then type:
```
export QBIN_DIR=$HOME/apps/Q6/bin   # absolute path to your local Q-binary folder (optional)
cd tests
./run_tests.sh
```
If `QBIN_DIR` is not set, pre-compiled Q binaries will be downloaded from Github.

#### Citations

To acknowledge the use of qtools in scientific publications, please specify the version and the
DOI of the latest release:  
[![DOI](https://zenodo.org/badge/80016679.svg)](https://zenodo.org/badge/latestdoi/80016679)  

*Example:*  
*...analysis was performed with qtools v0.5.10 (DOI: 10.5281/zenodo.842003).*

#### Bugs, feature requests, and contributions
Submit issues and feature requests to Github (https://github.com/mpurg/qtools/issues).  
To contribute create a pull request.

#### Author
Miha Purg (miha.purg@gmail.com)  

##### Contributors
Paul Bauer


