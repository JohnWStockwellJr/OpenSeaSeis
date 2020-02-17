#!/bin/bash

# -------------------------------------------------------
# SEASEIS make utility
# Windows 32bit version, compiled on Linux system
#
# NOTES
# This make utility works with the MinGW cross-compiler, installed on a Linux system.
# MinGW is available at http://mingw-w64.sourceforge.net/
# ...or through Linux distribution.
# For example, the Ubuntu packages for Windows cross-compilers are:
#   mingw-w64-x86-64-dev
#   mingw-w64-i686-dev
# (sudo apt-get install <package_name>)
#
# USAGE
# (1) Build Linux version using ./make_seaseis.sh
# (2) Build Windows 32bit version using ./make_seaseis_win32.sh
# (3) Copy directory 'win32_bin' to Windows system add to binary PATH
#
# -------------------------------------------------------
#
WINVER=32

# MinGW 'bin' directory must be in the binary path.

export SRCDIR=./src
export JAVADIR=./java
export WIN_LIBDIR=./win/win${WINVER}

#export CPP=i686-w64-mingw32-g++-win32
#export CC=i686-w64-mingw32-gcc-win32
#export LD=i686-w64-mingw32-g++-win32
#export F77=i686-w64-mingw32-gfortran-win32
export CPP=i686-w64-mingw32-g++
export CC=i686-w64-mingw32-gcc
export LD=i686-w64-mingw32-g++
export F77=i686-w64-mingw32-gfortran
export MAKE=make

export GLOBAL_FLAGS="-fexpensive-optimizations -O3 -Wno-long-long -Wall -pedantic"
export F77_FLAGS="-ffixed-line-length-132 -O3 -fexpensive-optimizations"
export COMMON_FLAGS="${GLOBAL_FLAGS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE"
export RM="rm -f"
export COPY=cp

export BUILD_FFTW=0 # Set to 1 if FFTW3 library is installed. Set details below where to find FFTW library
export BUILD_F77=1  # Set to 1 if Fortran compiler is available
export BUILD_SU=0   # Set to 1 if SU module shall be compiled. Requires special SU installation, see file README_SU
export BUILD_MPI=0  # Set to 1 to build with MPI enabled

export PLATFORM_WINDOWS=1

# -------------------------------------------------------------------

source $SRCDIR/make/linux/set_environment.sh

export BINDIR=${LIBDIR}/../win${WINVER}_bin
export OBJDIR=${LIBDIR}/../win${WINVER}_obj

echo "SeaSeis ${VERSION} source root directory:  '${CSEISDIR_SRCROOT}'"
echo "SeaSeis ${VERSION} obj/lib/bin root dir:   '${CSEISDIR}'"
echo "SeaSeis ${VERSION} library directory:      '${LIBDIR}'"

# -------------------------------------------------------------------
# Check if output directories exist or not. If not create them

if [ ! -d ${LIBDIR} ]; then
   echo "Directory ${LIBDIR} does not exist"
   echo "Run make_seaseis.sh first before running this make script"
   exit 1
fi

if [ ! -d $OBJDIR ]; then
  mkdir $OBJDIR
fi
if [ ! -d $OBJDIR ]; then
  echo "$OBJDIR is not a directory" ; exit 1
fi

if [ ! -d $BINDIR ]; then
  mkdir $BINDIR
fi
if [ ! -d $BINDIR ]; then
  echo "$BINDIR is not a directory" ; exit 1
fi

# -------------------------------------------------------------------
# Build CSEIS

# Create module 'h' file from module list (${SRCDIR}/include/cseis_modules.txt)
src/make/linux/prepare_cseis_build.sh

echo "Start building CSEIS..."

if [ -e $BINDIR/seaseis.exe ]; then
  ${RM} $BINDIR/seaseis.exe
fi

$MAKE -f $SRCDIR/make/win/Makefile_all
$MAKE -f $SRCDIR/make/win/Makefile_seaview

cp ${JAVADIR}/jar/CSeisLib.jar ${BINDIR}
cp ${JAVADIR}/jar/SeaView.jar ${BINDIR}
cp ${JAVADIR}/bin/seaview.bat ${BINDIR}
cp ${JAVADIR}/bin/plotimage.bat ${BINDIR}

echo "End..."
