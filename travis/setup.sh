#!/bin/sh

here=$(pwd)
fftw_dir="${here}/fftw"

echo "Compiling Charm++ (MPI)..."
wget http://charm.cs.illinois.edu/distrib/charm-6.10.1.tar.gz
tar xzf charm-6.10.1.tar.gz
mv charm-v6.10.1 charm-6.10.1
cd charm-6.10.1
env MPICXX=mpicxx ./build charm++ mpi-linux-x86_64 --with-production

echo "Testing Charm++ (MPI)..."
cd mpi-linux-x86_64/tests/charm++/megatest
make pgm
mpiexec -n 4 ./pgm

echo "Compiling FFTW..."
cd ../../../../..
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

echo "Done"