# quantum chem demo

## Setup

### OS
One of the central dependencies, `pyscf` is a python package that only works on MacOS and Linux. If you are on Windows, install [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) to use the Linux packages

Open a terminal for WSL or MacOS:

### Environment
Install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) if you do not already have it
To run the python files, create an environment with the necessary dependencies
`conda create -n stuy pyscf matplotlib rich matplotlib -c conda-forge`\
`conda activate stuy`\
`pip install pyberny`

### Terminal commands
Use the `cd` command to change directory to the one specified (e.g. `cd MetalHydroxides`)

## Running Python Simulations
You can run the python files as they are or explore some command-line flags available to adjust the simulation without needing to edit the code. If command-line flags are available, you can use the `-h` flag to list them. (e.g. `cd TitrationCurve` then `python titration.py -h`)


## Mathematica
For the kinetics simulation, get Wolfram Mathematica through your school's software store, if applicable. Otherwise it's expensive.
