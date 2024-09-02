#! /bin/bash

check_and_install_pkg() {
  REQUIRED_PKG=$1
  PKG_OK=$(dpkg-query -W --showformat='${Status}\n' $REQUIRED_PKG|grep "install ok installed")
  echo Checking for $REQUIRED_PKG: $PKG_OK
  if [ "" = "$PKG_OK" ]; then
    echo "No $REQUIRED_PKG. Setting up $REQUIRED_PKG."
    sudo apt-get --yes install $REQUIRED_PKG
  fi
}
sudo apt-get update
check_and_install_pkg "python3"

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
chmod a+x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create -n qreach python=3.9
conda activate qreach

export PYTHON_INCLUDE=`python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"`

# Install make
check_and_install_pkg "make"

# Download Boost C++ library
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz .
tar -xvf boost_1_81_0.tar.gz
cd boost_1_81_0/
export BOOST_PATH=$(pwd)
cd ..

# Install pip
check_and_install_pkg "python3-pip"

# Install invoke
pip install invoke

# Install pybind11
pip install pybind11

# Install gcc and g++
sudo apt-get update
check_and_install_pkg "build-essential"

# Install git
check_and_install_pkg "git"

# Install latex
check_and_install_pkg "texlive-latex-base"


# Clone repository
git clone git@github.com:Acdimy/qreach-tools.git
cd qreach-tools/

# Building QReach
cd python_pkg/
mkdir output/
# Not required; already fixed
# sed -i 's/-I${3}/-I{3}/g' tasks.py
invoke build-qreach
invoke build-pybind11

