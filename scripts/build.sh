#!/bin/sh

cd ..
root=$(pwd)
fftw_dir="${root}/fftw"

echo "Compiling Charm++ (MPI)..."
wget http://charm.cs.illinois.edu/distrib/charm-6.10.1.tar.gz
tar xzf charm-6.10.1.tar.gz
mv charm-v6.10.1 charm-6.10.1
cd charm-6.10.1
env MPICXX=mpicxx ./build charm++ mpi-linux-x86_64 --with-production

echo "Compiling FFTW..."
cd ..
wget http://www.fftw.org/fftw-2.1.5.tar.gz
tar xzf fftw-2.1.5.tar.gz
cd fftw-2.1.5
./configure --enable-float --enable-type-prefix --enable-static --prefix=${fftw_dir}
make
make install

echo "Downloading TCL libraries..."
cd ..
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64.tar.gz
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64-threaded.tar.gz
tar xzf tcl8.5.9-linux-x86_64.tar.gz
tar xzf tcl8.5.9-linux-x86_64-threaded.tar.gz
mv tcl8.5.9-linux-x86_64 tcl
mv tcl8.5.9-linux-x86_64-threaded tcl-threaded

echo "Compiling NAMD (MPI)..."
./config Linux-x86_64-g++ --charm-arch mpi-linux-x86_64
cd Linux-x86_64-g++
make -j $(nproc)

echo "Done"
