# Orthanq

Haplotype abundance quantification using a Bayesian model that is based on collecting evidences from Kallisto and Varlociraptor.

Orthanq is still a prototype and under active development.


## Building for development and manual installation

The easiest way to build Orthanq is via Conda/Mamba.
First, create the environment given by the `environment.yml` in this repository.

    mamba env create -n orthanq-dev -f environment.yml

Second, activate the environment with

    mamba activate orthanq-dev

Third, build Orthanq via

    HDF5_DIR=$CONDA_PREFIX RUSTFLAGS='-L $CONDA_PREFIX/include' cargo build

    
