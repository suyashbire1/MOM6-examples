#!/bin/bash -f

source $MODULESHOME/init/bash
#module load anaconda
#module load nco/4.5.4
#module load mpich2/1.2.1p1
#module load intel_compilers

source activate MIDAS

echo -n Started $0 in ; pwd


# Installing this file avoids a slow and unreliable file server (that would otherwise be ftp'd from ftp.oce.orst.edu)
cp -n /archive/gold/datasets/obs/tpxo7_atlas_netcdf.tar.Z .

# Work around for environment problem inside MIDAS
#setenv PYTHONPATH $cwd/local/lib

# Run through the work flow
make
