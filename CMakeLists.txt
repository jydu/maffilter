# CMake script for GenomeTools package
# Author: Julien Dutheil
# Created: 27/04/2010

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(testnh C CXX)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF()

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF()

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})

  IF(BUILD_STATIC)
    FIND_PATH(${lib_name}_LIBRARY_DIR lib${lib_name}.a HINTS ${CMAKE_PREFIX_PATH} ${CMAKE_PREFIX_PATH}/lib)
    IF(${lib_name}_LIBRARY_DIR)
      SET(${lib_name}_LIBRARY ${${lib_name}_LIBRARY_DIR}/lib${lib_name}.a)
    ELSE(${lib_name}_LIBRARY_DIR)
      MESSAGE(FATAL_ERROR "${lib_name} static library could not be found.")
    ENDIF(${lib_name}_LIBRARY_DIR)
  ELSE(BUILD_STATIC)
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  ENDIF(BUILD_STATIC)

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   libraries: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

#Find the libraries. The order is very important for staitc linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Seq/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-numcalc NumCalc/NumTools.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-utils Utils/Clonable.h)

# Find the zlib installation
FIND_PACKAGE(ZLIB REQUIRED)
INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
SET(LIBS ${LIBS} ${ZLIB_LIBRARIES})

# Find the bz2 installation
FIND_PACKAGE(BZip2 REQUIRED)
INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
SET(LIBS ${LIBS} ${BZIP2_LIBRARIES})

#Find boost libraries
SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.34.0 COMPONENTS iostreams )

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

# Subdirectories
ADD_SUBDIRECTORY(GenomeTools)

# Packager
SET(CPACK_PACKAGE_NAME "genometools")
SET(CPACK_PACKAGE_VENDOR "Julien Dutheil")
SET(CPACK_PACKAGE_VERSION "0.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Tools for genomics")
SET(CPACK_RESOURCE_FILE_LICENSE "COPYING")
SET(CPACK_RESOURCE_FILE_AUTHORS "AUTHORS")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.svn"
 ".*\\\\.gz"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 "\\\\..*\\\\.swp"
 "GenomeTools/\\\\..*"
 "GenomeTools/mafinfo"
 "debian/tmp"
 "debian/mafinfo"
 ${CPACK_SOURCE_IGNORE_FILES}
)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(deb fakeroot debian/rules binary)
ADD_DEPENDENCIES(deb dist info)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info)

ENDIF()

