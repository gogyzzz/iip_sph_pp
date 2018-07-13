#!/bin/sh

#Remove existing cmake
echo "Remove Existing Cmake"
sudo apt remove cmake
sudo apt purge --auto-remove cmake

echo "Downloading cmake 3.10.3"
wget https://cmake.org/files/v3.10/cmake-3.10.3.tar.gz
tar -xzvf cmake-3.10.3.tar.gz
cd cmake-3.10.3
./bootstrap
make
sudo make install
sudo cp bin/cmake /usr/bin  
sudo cp bin/cpack /usr/bin  
sudo cp bin/ctest /usr/bin  

echo "CMAKE INSTALL COMPLETE"
cmake --version
