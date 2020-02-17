#!/bin/bash

# Usage:   make_seaseis.sh [option]
# Options: clean, bleach, verbose, debug

#*****************************************************
# Compiler settings
#
export BUILD_FFTW=0 # Set to 1 if FFTW3 library is installed. Set details below where to find FFTW library
export BUILD_F77=1  # Set to 1 if Fortran compiler is available
export BUILD_SU=0   # Set to 1 if SU module shall be compiled. Requires special SU installation, see file README_SU
export BUILD_MPI=0  # Set to 1 to build with MPI enabled

if [ ${BUILD_MPI} -eq 1 ]; then
    export CPP=mpiCC
else
    export CPP=g++
fi
export F77=gfortran
export LD=$CPP

#*****************************************************
# Make file settings, command line arguments
#
export MAKE_BLEACH=0
export MAKE_CLEAN=0
export MAKE_DEBUG=0
export VERBOSE=0
export NO_MAKE_BUILD=0

# Determine special setting for shared object build (SONAME)
unameString=$( uname | awk '{print tolower($0)}' )
if [ $unameString == "linux" ]; then
  export platform="linux"
  export SONAME=soname
elif [ $unameString == "darwin" ]; then
  export platform="apple"
  export SONAME=install_name
else
  # ..otherwise just assume Linux
  export platform="linux"
  export SONAME=soname
fi

make_argument=""
# Uncomment to make 'silent'
# make_argument="-s"

export GLOBAL_FLAGS="-fexpensive-optimizations -O3 -Wno-long-long -Wall -pedantic -Wformat-overflow=0"
export F77_FLAGS="-ffixed-line-length-132"
export RM="rm -f"

for arg in $@
do
    if [ $arg == "verbose" ]; then
        export VERBOSE=1
        make_argument=""
    fi
    if [ $arg == "debug" ]; then
        export MAKE_DEBUG=1
        make_argument=""
        GLOBAL_FLAGS="-g -DOS_DEBUG=1 -Wno-long-long -Wall -pedantic -Wconversion"
        F77_FLAGS+=" -g"
    fi
    if [ $arg == "clean" ]; then
        export MAKE_CLEAN=1
        make_argument="clean"
    fi
    if [ $arg == "bleach" ]; then
        export MAKE_BLEACH=1
        make_argument="bleach"
    fi
done

numMakeOptions=$(echo "${MAKE_BLEACH} + ${MAKE_CLEAN} + ${MAKE_DEBUG}" | bc -l)
if [ ${numMakeOptions} -gt 1 ]; then
    echo "ERROR: Too many make options. Specify one option only: debug, clean or bleach"
    exit
fi

#********************************************************************************
# Setup environment variables for Seaseis directories
#

source src/make/linux/set_environment.sh

# *** FOR DEVELOPERS:  Comment out the following three lines to skip legal statement etc
./license.sh
#./mailhome.sh
export VERBOSE=1

#********************************************************************************
# Check if FFTW3 library is installed
# How to install fftw library for use with OpenSeaSeis:
# 1) Download and extract fftw source code distribution from http://www.fftw.org/download.html
# 2) ./configure --enable-float --enable-shared 
# 3) make
# 4) make install

export INCDIR_FFTW=/usr/local/include
export LIBDIR_FFTW=/usr/local/lib
export LIBFFTW=fftw3f

if [ ${BUILD_FFTW} -eq 1 ]; then
  export FOUND_FFTW=$($CPP -l${LIBFFTW} -L${LIBDIR_FFTW} 2>&1 | grep -i "cannot find -l" | wc -l | awk '{print "1-"$1}' | bc )
  if [ ${FOUND_FFTW} -eq 0 ]; then
      echo "WARNING: Library 'lib${LIBFFTW}.so' not found in library path."
      echo " - Seaseis FFT modules will be compiled without FFTW."
      echo " - Visit www.fftw.org to download and install FFTW libraries."
      export BUILD_FFTW=0
  fi
fi

echo "SeaSeis ${VERSION} source root directory:  '${CSEISDIR_SRCROOT}'"
echo "SeaSeis ${VERSION} obj/lib/bin root dir:   '${CSEISDIR}'"
echo "SeaSeis ${VERSION} library directory:      '${LIBDIR}'"
if [ ${BUILD_FFTW} -eq 1 ]; then
  echo "FFTW library:      '${LIBDIR_FFTW}/lib${LIBFFTW}.so'"
fi

#********************************************************************************
# Make Seaseis
#

export COMMON_FLAGS="${GLOBAL_FLAGS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE"
export LD_FLAGS=${GLOBAL_FLAGS}
export LIB_FLAGS="-lc"
if [ ${BUILD_FFTW} -eq 1 ]; then
   export COMMON_FLAGS="${COMMON_FLAGS} -D USE_FFTW -I${INCDIR_FFTW}"
   export LIB_FLAGS="-Wl,-rpath,${LIBDIR_FFTW} -L${LIBDIR_FFTW} -lc -l${LIBFFTW}"
fi
if [ ${BUILD_MPI} -eq 1 ]; then
   export COMMON_FLAGS="${COMMON_FLAGS} -D USE_MPI"
fi

${CSEISDIR_SRCROOT}/src/make/linux/cmake.sh ${make_argument}

if [ ${MAKE_BLEACH} -eq 0 -a ${MAKE_CLEAN} -eq 0 ]; then
    echo "Auto-generate HTML self-documentation:  ${BINDIR}/seaseis -html > ${DOCDIR}/SeaSeis_help.html"
    ${BINDIR}/seaseis -html > ${DOCDIR}/SeaSeis_help.html
    ${THISDIR}/set_links ${VERSION}
fi

echo "Done."
