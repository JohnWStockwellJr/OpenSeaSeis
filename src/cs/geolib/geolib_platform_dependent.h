/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CSEIS_PLATFORM_DEPENDENT_H
#define CSEIS_PLATFORM_DEPENDENT_H

/**
* This is the 'only' platform dependent file in CSeis.
* Include in files where distinctions need to be made between different platforms.
* For SGI and Solaris, the following define statements have not been tested yet
*/


#if defined(__linux__) || defined(__linux)
 #define PLATFORM_LINUX 1

#elif defined(__APPLE__) || defined(__APPLE)
 #define PLATFORM_APPLE 1
 #define off64_t off_t
 #define fseeko64 fseeko
 #define ftello64 ftello
 
#elif defined(__sun__) || defined(__sun)
 #define PLATFORM_SOLARIS 1

#elif defined(__sgi__) || defined(__sgi)
 #define PLATFORM_SGI 1

#elif defined(WIN32)
 #define PLATFORM_WINDOWS 1

#endif

#if defined(__itanium__)
 #define ARCHITECTURE_ITANIUM 1
#endif

/*
 * Apparently, the Gnu g++ compiler works with either separator, on both Windows and Unix systems
#ifndef PLATFORM_WINDOWS
static char const DIR_SEPARATOR = '\\';
#else
static char const DIR_SEPARATOR = '/';
#endif
*/

#endif

