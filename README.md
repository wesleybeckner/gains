
![](gains.png)

## Genetic Algorithm for Identifying Novel Structures
[![Build Status](https://travis-ci.org/wesleybeckner/gains.svg?branch=master)](https://travis-ci.org/wesleybeckner/gains)
[![PyPI version](https://badge.fury.io/py/gains.svg)](https://badge.fury.io/py/gains)
[![Coverage Status](https://coveralls.io/repos/github/wesleybeckner/gains/badge.svg?branch=master)](https://coveralls.io/github/wesleybeckner/gains?branch=master)
========

GAINS — Genetic Algorithm for Identifying Novel Structures — is a project
that enables molecular design and computational screening of
small molecules. Built on the molecular functionality of RDKit, GAINS is employable across a spectrum of small-molecule design problems.

## Installation

### Dependencies

GAINS requires:

* python (>= 3.6)
* scikit-learn (>= 0.19.1)
* rdkit (>= 2017.09.1)
* salty-ilthermo (>= 0.2)

Note that scikit-learn 0.18.1 will raise a warning when loading in property models to the engine.

To take full advantage of rdkit you will also need Matplotlib >= 1.3.1.

#### User installation

You will first need to install [rdkit](http://www.rdkit.org/docs/GettingStartedInPython.html):
```
conda create -n py36 python=3.6 anaconda
# activate the new virtual environment, e.g. on OSX/Linux
source activate py36
# on Windows
# activate py36
conda install -c rdkit rdkit
```

GAINS can then be installed with:
```
pip install gains
```

## Development

GAINS is currently underdevelopment by researchers at the University of Washington. Our  research page can be found [here](http://www.prg.washington.edu).

### Testing

After installation, you can launch the test suite from outside the source directory (you will need to have the pytest package installed):
```
pytest gains
```
