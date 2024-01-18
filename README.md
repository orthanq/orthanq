[![DOI](https://zenodo.org/badge/475406908.svg)](https://zenodo.org/badge/latestdoi/475406908)

# ![Orthanq](logo.svg)

Haplotype abundance quantification using a Bayesian model that is based on collecting evidences from Varlociraptor.

Orthanq is under active development.

## Building for development and manual installation

The easiest way to build Orthanq is via Conda/Mamba.
First, create the environment given by the `environment.yml` in this repository.

    mamba env create -n orthanq-dev -f environment.yml

Second, activate the environment with

    mamba activate orthanq-dev

Third, build Orthanq via

    HDF5_DIR=$CONDA_PREFIX RUSTFLAGS='-L $CONDA_PREFIX/include' cargo build

Last, conda lib location is set via

    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    
