mkdir build                                                                                                                                                                                                                                                                                                                                                                                                               
cd build
source $DD4hep_DIR/bin/thisdd4hep.sh

cmake .. -DCMAKE_INSTALL_PREFIX=../InstallArea -DBoost_NO_BOOST_CMAKE=ON
make -j4 install
cd ../InstallArea
# source bin/thislcgeo.sh
source bin/thisk4geo.sh

export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib:$PYTHONPATH
cd ..
export k4geo_DIR=$(pwd)

export K4GEO=$(pwd)
