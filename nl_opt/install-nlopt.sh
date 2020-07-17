~!/bin/bash

# Usage:
# bash install-nlopt.sh download_and_install_directory_for_nlopt

DOWNLOAD_DIR=$1
wget -P ${DOWNLOAD_DIR} https://github.com/stevengj/nlopt/archive/v2.6.2.tar.gz
cd ${DOWNLOAD_DIR}
tar -xzvf v2.6.2.tar.gz
cd ${DOWNLOAD_DIR}/nlopt-2.6.2
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${DOWNLOAD_DIR}/nlopt ..
make
make install
