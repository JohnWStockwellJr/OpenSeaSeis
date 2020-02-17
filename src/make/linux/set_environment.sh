#!/bin/bash

#********************************************************************************
# Setup environment variables for Seaseis directories
#
export THISDIR=$(pwd)
export CSEISDIR_SRCROOT=${THISDIR}
export CSEISDIR_LIBROOT=${THISDIR}/..

# Retrieve version number
file_ver_chk=${CSEISDIR_SRCROOT}/src/cs/system/cseis_defines.h
if [ ! -e ${file_ver_chk} ]; then
  echo "File not found: ${file_ver_chk}"
  exit -1
fi
export VERSION=$(grep CSEIS_VERSION ${file_ver_chk} | awk 'BEGIN{FS="\""}{print "v"$2}')
export CSEISDIR=${CSEISDIR_LIBROOT}/lib_${VERSION}
export SRCDIR=${CSEISDIR_SRCROOT}/src
export JAVADIR=${CSEISDIR_SRCROOT}/java
export LIBDIR=${CSEISDIR}/lib
export BINDIR=${CSEISDIR}/bin
export OBJDIR=${CSEISDIR}/obj
export DOCDIR=${CSEISDIR}/doc
export INCDIR=${LIBDIR}/include
mkdir -p ${CSEISDIR}
mkdir -p ${INCDIR}
