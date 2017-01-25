### README 

#### Description

A Python interface library (Qpyl) and a bunch of useful command line tools for [Q](http://xray.bmc.uu.se/~aqwww/q/), your friendly EVB simulation package.

###### Useful?

- Library/parameter file conversion (oplsaa, amber).
- FEP file generation for EVB simulations.
- Input generation (MD, FEP, EVB).
- EVB mapping and automatic fitting of EVB parameters.
- Output parsing/analysis (MD, FEP, EVB).

###### Tests?

Working on it, jeez...

###### Documentation?

*Please* contribute, I pay in beer/pizza/karma points...


#### Requirements

The tools pride themselves on not having any dependencies, except *matplotlib*
which is required only for *q_plot.pt*, a non-essential utility used mostly for
debugging.

The only requirement is *Python, version 2.7*

*Note: The tools will likely fail misserably on Windows,
since development (and testing) is done exclusively on GNU/Linux.*

#### Installation (\*nix, Bash)
Clone this repository to your favourite directory:  
```
mkdir -p ~/bin && cd ~/bin
git clone https://bitbucket.org/grupim/qscripts.git
```



Add this to your `~/.bashrc` or `~/.bash_profile`:  
```
source $HOME/bin/qscripts/qscripts_init.sh
```

Run the config script:  
```
qscripts_config.py
```



#### Other
Bugs, suggestions and questions: miha.purg@gmail.com


###### Contributors
Paul Bauer 



