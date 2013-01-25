#!/bin/bash

if [ -z "$1" ] ; then
    echo "usage: $0 <src-tarball>"
    exit 1
fi

CMAKE=/home/lwilcox/opt/cmake/2.8.10.2/bin/cmake

TARBALL_BASE=`basename $1`
TARBALL_ABSPATH=$( cd $(dirname $1); pwd)/$(basename $1)
BFAM_VERSION=`echo $TARBALL_BASE | sed 's/bfam-\(.*\)-Source.tar.gz/\1/'`
BUILD_DIR=$HOME/work/builds/bfam/$BFAM_VERSION
SRC_BASE=$HOME/work/src
SRC_DIR=$SRC_BASE/bfam-$BFAM_VERSION-Source

echo "Make: $SRC_BASE"
mkdir -p $SRC_BASE

pushd $SRC_BASE
echo "Extracting: $TARBALL_ABSPATH"
tar -xzvf $TARBALL_ABSPATH
popd

module unload compile mpi

module load compile/intel/13.0 mpi/openmpi/1.6.3
FLAVOR=intel-openmpi-debug
echo "Building: $FLAVOR"
rm -rf $BUILD_DIR/$FLAVOR
mkdir -p $BUILD_DIR/$FLAVOR
pushd $BUILD_DIR/$FLAVOR
$CMAKE -DCMAKE_BUILD_TYPE="Debug" $SRC_DIR
make -j 8
make test
popd
module unload compile mpi

module load compile/intel/13.0 mpi/openmpi/1.6.3
FLAVOR=intel-openmpi-release
echo "Building: $FLAVOR"
rm -rf $BUILD_DIR/$FLAVOR
mkdir -p $BUILD_DIR/$FLAVOR
pushd $BUILD_DIR/$FLAVOR
$CMAKE -DCMAKE_BUILD_TYPE="Release" $SRC_DIR
make -j 8
make test
popd
module unload compile mpi

module load compile/intel/13.0 mpi/mvapich2/1.9a2
FLAVOR=intel-mvapich2-debug
echo "Building: $FLAVOR"
rm -rf $BUILD_DIR/$FLAVOR
mkdir -p $BUILD_DIR/$FLAVOR
pushd $BUILD_DIR/$FLAVOR
$CMAKE -DCMAKE_BUILD_TYPE="Debug" $SRC_DIR
make -j 8
make test
popd
module unload compile mpi

module load compile/intel/13.0 mpi/mvapich2/1.9a2
FLAVOR=intel-mvapich2-release
echo "Building: $FLAVOR"
rm -rf $BUILD_DIR/$FLAVOR
mkdir -p $BUILD_DIR/$FLAVOR
pushd $BUILD_DIR/$FLAVOR
$CMAKE -DCMAKE_BUILD_TYPE="Release" $SRC_DIR
make -j 8
make test
popd
module unload compile mpi

## module load compile/pgi/12.10 mpi/openmpi/1.6.3
## FLAVOR=pgi-openmpi-debug
## echo "Building: $FLAVOR"
## rm -rf $BUILD_DIR/$FLAVOR
## mkdir -p $BUILD_DIR/$FLAVOR
## pushd $BUILD_DIR/$FLAVOR
## $CMAKE -DCMAKE_BUILD_TYPE="Debug" $SRC_DIR
## make -j 8
## make test
## popd
## module unload compile mpi
##
## module load compile/pgi/12.10 mpi/openmpi/1.6.3
## FLAVOR=pgi-openmpi-release
## echo "Building: $FLAVOR"
## rm -rf $BUILD_DIR/$FLAVOR
## mkdir -p $BUILD_DIR/$FLAVOR
## pushd $BUILD_DIR/$FLAVOR
## $CMAKE -DCMAKE_BUILD_TYPE="Release" $SRC_DIR
## make -j 8
## make test
## popd
## module unload compile mpi

echo "Done building"
