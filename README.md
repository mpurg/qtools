[![Travis](https://api.travis-ci.org/mpurg/qtools.svg?branch=master)](https://travis-ci.org/mpurg/qtools)  [![Coverage Status](https://coveralls.io/repos/github/mpurg/qtools/badge.svg?branch=master)](https://coveralls.io/github/mpurg/qtools?branch=master)

  
  

### README 

#### Description

qtools is a toolset comprised of a Python library Qpyl and a set of command-line tools, for use with the [molecular dynamics simulation package *Q*](http://xray.bmc.uu.se/~aqwww/q/). This open-source project was created by Miha Purg (currently at Uppsala University) in order to provide Q-users with a set of well-tested tools for automating repetitive and/or tedious tasks, debugging simulations, and analysing results. The author strives to achieve a good balance between the level of user-control over *Q* (i.e., relatively low abstraction level) and the ease-of-use provided by qtools.

*Disclaimer: qtools is not affiliated with *Q* nor should it be considered an official toolset for Q.*


#### Functionality

- Python API for interfacing with Q executables.
- EVB mapping and automatic fitting of off-diagonal and gas shift parameters.
- Support for quantum classical path and group exclusion calculations.
- Linear response approximation (LRA) calculations (including group contributions)
- Parameter conversions (OPLS, Amber/GAFF).
- FEP file generation for EVB simulations.
- Input generation for MD/FEP simulations and EVB calculations.
- Output parsing and analysis of MD, FEP and EVB.
- Graphical tool for visualizing results and debugging simulations:
  - Free energy profiles
  - LRA, REORG, group contributions
  - Energy breakdown (VDW, El, Bonds, ...)


#### Requirements

The library and tools do not have any dependencies, except *matplotlib*
which is required only for *q_plot.py*, a non-essential utility used mostly for
troubleshooting simulations.

The only requirement is *Python, version 2.7*

*Note: The command-line tools were designed to work in a typical high-perfomance-computer Unix environment, and have not been tested with Microsoft Windows.*


#### Installation (Linux)

Clone the repository to your local directory:  
```
mkdir -p ~/bin && cd ~/bin
git clone https://github.com/mpurg/qtools
```

Add this line to your `~/.bashrc` or `~/.bash_profile`:  
```
source $HOME/bin/qtools/qtools_init.sh
```

Run the CLI config script:  
```
qscripts_config.py
```

#### Documentation

The starting point should be the `docs/tutorials` folder.
It contains working examples of how you can use the command line tools to prepare, simulate and analyse a simple SN2 reaction with EVB. Additionally, a recently published chapter in the serial Methods in Enzymology provides step-by-step instructions for performing EVB simulations on a real enzymatic system with `qtools`. Manuscript available [here](https://doi.org/10.1016/bs.mie.2018.06.007) (paywalled).
  
For those inclined on using the Python API directly, jupyter-notebook examples calling the Qpyl API can be found in `docs/examples`. The API reference can be found in `docs/api/html`. Additionally, the tests, and the command-line interface (qscripts-cli) can be useful to explore the API.
  
#### Testing

Automated Testing is performed on [Travis-CI](https://travis-ci.org/mpurg/qtools),
with [pytest](https://docs.pytest.org/en/latest/) for Qpyl
([code coverage](https://coveralls.io/github/mpurg/qtools?branch=master)).
and simple regression tests for CLI tools.
To run the tests locally, make sure you have `pytest` installed and simply type:
```
cd tests
./run_tests.sh
```

#### Citations

To acknowledge the use of qtools in scientific publications, please specify the version and the
DOI of the latest release:  
[![DOI](https://zenodo.org/badge/80016679.svg)](https://zenodo.org/badge/latestdoi/80016679)  

*Example:*  
*...analysis was performed with qtools v0.5.10 (DOI: 10.5281/zenodo.842003).*


#### Author
Miha Purg (miha.purg@gmail.com)  

##### Contributors
Paul Bauer

