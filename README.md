[![Travis](https://api.travis-ci.org/mpurg/qtools.svg?branch=master)](https://travis-ci.org/mpurg/qtools)  [![Coverage Status](https://coveralls.io/repos/github/mpurg/qtools/badge.svg?branch=master)](https://coveralls.io/github/mpurg/qtools?branch=master)

  
  

### README 

#### Description

qtools is a Python interface library (Qpyl) and a set of useful
command line tools for working with the [*Q* molecular dynamics 
simulation package](http://xray.bmc.uu.se/~aqwww/q/).
This open-source project was created by Miha Purg
(currently at Uppsala University) in order to provide Q-users with a set of well tested
tools for automating repetitive and/or tedious tasks, debugging simulations, and analysing results.
The author strives to achieve a good balance between the level of
user-control over *Q* (i.e., relatively low abstraction level) and the ease-of-use provided by qtools.

*Disclaimer: qtools is not affiliated with *Q* nor should it be considered an official toolset for Q.*


#### Available features

- Converting parameters to Q format (OPLS, Amber/GAFF).
- FEP file generation for EVB simulations.
- Input generation for MD/FEP simulations and EVB calculations.
- EVB mapping and automatic fitting of off-diagonal and gas shift parameters.
- Output parsing and analysis of MD, FEP and EVB.
- Graphical tool for visualizing results and debugging simulations:
  - Free energy profiles
  - LRA, REORG, Group contributions
  - Energy breakdown (VDW, El, Bonds, ...)


#### Requirements

The tools pride themselves on not having any dependencies, except *matplotlib*
which is required only for *q_plot.py*, a non-essential utility used mostly for
debugging.

The only requirement is *Python, version 2.7*

*Note: The tools will likely fail misserably on Windows,
since development (and testing) is done exclusively on GNU/Linux.*
  



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

Run the config script:  
```
qscripts_config.py
```

#### Documentation

The starting point should be the `docs/tutorials` folder.
For now, it contains one working example of how you can use the command line tools to prepare, simulate and analyse a simple SN2 reaction with EVB.  
  
For python developers, jupyter-notebook examples of interfacing Qpyl can be found in `docs/examples`. Additionally the tests, qscripts_cli and the Sphinx API reference (in `docs/manual`) can be useful to explore the API.

More documentation will be added as soon as possible.


#### Testing

Testing is done [here](https://travis-ci.org/mpurg/qtools),
with [pytest](https://docs.pytest.org/en/latest/) for Qpyl
([code coverage](https://coveralls.io/github/mpurg/qtools?branch=master)).
and simple regression tests for CLI tools.
To run the tests locally, make sure you have `pytest` installed and simply type:
```
cd tests
./run_tests.sh
```

*Note:
Parts of the code rely on (at the moment still proprietary) Q
binaries and are thus tested elsewhere.*


#### Citations

To acknowledge the use of qtools in scientific publications, please specify the version and the
DOI of the latest release:  
[![DOI](https://zenodo.org/badge/80016679.svg)](https://zenodo.org/badge/latestdoi/80016679)  

*Example:*  
*...analysis was performed with qtools v0.5.10 (DOI: 10.5281/zenodo.842003).*


#### Author
Miha Purg (miha.purg@gmail.com)  

##### Contributors
Paul Bauer (paul.bauer@icm.uu.se)  

