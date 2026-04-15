Many-body Dipole and Polarizability codes from MB-MD
====================================================

On the master branch, you can find programs that calculate the total dipole and
polarizability in terms of a many-body expansion. As described in
[this paper](http://pubs.acs.org/doi/abs/10.1021/ct501131j), this is implemented as a sum of
one-body (1B) properties, short-ranged 2B properties described (in part) through
permutationally invariant polynomials, and long-ranged 2B plus higher-order
many-body contributions coming from an inducible dipole electrostatics model.

For calculating the dipole/polarizability, two types of programs can be found
to process

- Bulk simulation trajectories from dlpoly `POSITION_CMD` files (code prefix: `bulk_property_model.cpp`)
- Cluster configurations in xyz format (code prefix: `cluster_property_model.cpp`)

As is hopefully suggested by this notation, the format of the file name is `X_Y_Z.cpp`, where `X` indicates the type of configuration to be processed
(bulk or cluster), `Y` indicates whether the dipole or polarizability will be
calculated, and `Z` indicates the model that is used to calculate the property
(e.g., q-TIP4P/f or MB-$\mu$ for the dipole moment).

# Workflow for calculating vibrational spectra

Other codes can also be found in this repo, notably, programs for the
calculation of vibrational spectra. Vibrational spectra can be obtained through
three steps: 1) calculating the time dependence of an electric property
(e.g., the dipole moment), 2) calculating the (ensemble-averaged) time
correlation function of these properties (e.g., the time autocorrelation
autocorrelation function of the dipole moment, for the IR spectra), and 3)
calculating the Fourier transform of the time correlation function to give
the vibrational spectra. In this repository, each of these three steps is
implemented in its own program.

  1. `bulk_dip` calculates the time dependence of the dipole moment using the simulation trajectory from DLPOLY.
  2. `response_ir` calculates the time autocorrelation function (response function) of the dipole moment using the output of step 1.
  3. `spectrum-ir` calculates the infrared spectrum using the output of step 2.

# Polynomial fitting for dipole/polarizability

These codes are contained in a separate branch of this repository, called `clean_fit`

# Compilation directions
to compile, run the following:

```
./autogen.sh
./configure CXX=mpicxx CC=mpicc --enable-maintainer-mode --disable-silent-rules PKG_CONFIG_PATH=/usr/lib64/pkgconfig CFLAGS='-g -fopenmp -O0 -Wall -std=c99' CXXFLAGS='-g -fopenmp -O0 -Wall'
./config.status
make -j4
```
# mb-spec-dev-dev
