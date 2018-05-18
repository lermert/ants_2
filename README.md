# ants_2
Download, process and correlate continuous seismic data 

# Documentation
Information on usage can be found in the ants-doc.tar archive. After installation (see below), start with ants_quickstart.pdf.

# Dependencies
obspy

mpi4py (uses openmpi or mpich)

h5py (uses hdf5)

basemap

mock

click

pyasdf (https://github.com/SeismicData/pyasdf)

It is convenient to install the dependencies with anaconda (https://www.continuum.io/downloads). To install obspy and pyasdf, adding the conda-forge channel is necessary, see  https://github.com/obspy/obspy/wiki/Installation-via-Anaconda. If you want to write correlations to pyasdf directly, it is currently necessary to install it with parallel hdf5 support; see http://seismicdata.github.io/pyasdf/installation.html. In most cases it is easier to write them to SAC and convert to pyasdf after all correlations are computed.


# Installation

After cloning, ants should be set up by running:

pip install -v -e .

in the ants_2 directory.
If successfully installed, one should be able to call the scripts with:

ants

This should display a list of commands.

# mpi4py

Sometimes mpi4py causes problems when running on mac OS. In some cases installing mpi4py with pip after having installed the other dependencies with conda works better.
