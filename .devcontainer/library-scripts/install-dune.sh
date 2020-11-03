#!/usr/bin/env bash

# requirements
apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y zlib1g zlib1g-dev paraview gmsh  bison flex libmetis-dev libparmetis-dev \
    libalberta-dev libsuperlu-dev libarpack++2-dev libsuitesparse-dev libpsurface-dev libgmp-dev libhypre-dev
# texlive-base texlive-latex-recommended texlive texlive-latex-extra \
#    texlive-bibtex-extra texlive-fonts-extra texlive-full doxygen 
    
sed -i 's/rights="none" pattern="EPS"/rights="read|write" pattern="EPS"/g' /etc/ImageMagick-6/policy.xml

# install dune
mkdir -p /tmp/dune
cd /tmp/dune || exit 1

VERSION=2.6.0

# core modules
wget https://dune-project.org/download/"${VERSION}"/dune-common-"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-geometry-"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-grid-"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-grid-howto-"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-istl-"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-localfunctions-"${VERSION}".tar.gz

tar -xvzf dune-common-"${VERSION}".tar.gz
tar -xzf dune-geometry-"${VERSION}".tar.gz
tar -xzf dune-grid-"${VERSION}".tar.gz
tar -xzf dune-grid-howto-"${VERSION}".tar.gz
tar -xzf dune-istl-"${VERSION}".tar.gz
tar -xzf dune-localfunctions-"${VERSION}".tar.gz

# extension modules
#wget https://gitlab.dune-project.org/extensions/dune-alugrid/-/archive/v"${VERSION}"/dune-alugrid-v"${VERSION}".tar.gz
#tar -xzf dune-alugrid-v"${VERSION}".tar.gz
wget https://dune-project.org/download/"${VERSION}"/dune-uggrid-"${VERSION}".tar.gz
tar -xzf dune-uggrid-"${VERSION}".tar.gz
#git clone https://git.imp.fu-berlin.de/agnumpde/dune-solvers.git
#git clone https://gitlab.dune-project.org/extensions/dune-alugrid

dune-common-"${VERSION}"/bin/dunecontrol configure -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_DISABLE_FIND_PACKAGE_MPI=1
dune-common-"${VERSION}"/bin/dunecontrol make install


# rm -r /tmp/dune
