#!/usr/bin/env bash

set -e

mkdir -p /tmp/deps
cd /tmp/deps || exit 1

# ITSOL
wget http://www-users.cs.umn.edu/~saad/software/ITSOL/ITSOL_2.tar.gz
tar xzf ITSOL_2.tar.gz

cd ITSOL_2 || exit 1
make lib
cp INC/* /usr/local/include
cp LIB/* /usr/local/lib
cd -
rm ITSOL_2.tar.gz
rm -r ITSOL_2

# HYPRE
wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.11.2.tar.gz
tar xzf hypre-2.11.2.tar.gz
cd hypre-2.11.2/src/cmbuild || exit 1
cmake .. -GNinja -DCMAKE_INSTALL_PREFIX=/usr/local
ninja install
cd -
rm hypre-2.11.2.tar.gz
rm -r hypre-2.11.2

