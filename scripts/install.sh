### Source script to install k4geo 

# using the nightly release
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh                                                                                                                                                                                                                                                                                                              

script_folder="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Get the location of this script, to execute correctly other commands
cd $script_folder/..

echo "Running on xml file $xml"
mkdir build
cd build
source $DD4hep_DIR/bin/thisdd4hep.sh

cmake .. -DCMAKE_INSTALL_PREFIX=../InstallArea -DBoost_NO_BOOST_CMAKE=ON
make -j4 install
cd ../InstallArea
source bin/thislcgeo.sh

export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib:$PYTHONPATH
cd ..
export k4geo_DIR=$(pwd)