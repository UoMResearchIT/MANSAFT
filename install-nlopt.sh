cd ${HOME}/Downloads
wget https://github.com/stevengj/nlopt/archive/v2.6.2.tar.gz
tar -xzvf v2.6.2.tar.gz
cd nlopt-2.6.2
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/Downloads/nlopt ..
make
make install
