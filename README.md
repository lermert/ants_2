# ants_2
Download, process and correlate continuous data 

# Info on usage
Information on usage can be found in the ants-doc.tar archive. After installation (see below), start with ants_quickstart.pdf.

# Dependencies
obspy

mpi4py (uses openmpi or mpich)

h5py (uses hdf5)

basemap

mock

click

It is recommended to install the dependencies with anaconda (https://www.continuum.io/downloads), except in the case of mpi4py where depending on OS installation using pip may work better.

# Installation

After cloning, ants should be set up by running:

pip install -v -e .

in the ants_2 directory.
If successfully installed, one should be able to call the scripts with:

ants

This should display a list of commands.

# mpi4py

Sometimes mpi4py causes problems when running on mac OS. In some cases installing mpi4py with pip after having installed the other dependencies with conda works better.
