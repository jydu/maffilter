# CMake script for the MafFilter package
# Author: Julien Dutheil
# Created: 27/04/2010

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(maffilter C CXX)

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

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NOT DOC_COMPRESS)
  SET(DOC_COMPRESS gzip CACHE STRING
      "Set program for compressing documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS)

IF(NOT DOC_COMPRESS_EXT)
  SET(DOC_COMPRESS_EXT gz CACHE STRING
      "Set extension of compressed documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS_EXT)


IF(NO_DEP_CHECK)
  MESSAGE("-- Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc -static-libstdc++ ${CMAKE_CXX_FLAGS}")
ENDIF()

#build info?
IF(NOT DEFINED INFO)
  SET(INFO UNIX)
ENDIF(NOT DEFINED INFO)

#build man pages?
IF(NOT DEFINED MAN)
  SET(MAN UNIX)
ENDIF(NOT DEFINED MAN)

#find executables for documentation
IF(MAN)
  FIND_PROGRAM(NROFF_EXE NAMES nroff)
  IF(NROFF_EXE)
    MESSAGE("-- Found nroff here: ${NROFF_EXE}")
    MESSAGE("   Adding targets: man")

    ADD_CUSTOM_TARGET(man
      ALL
      COMMAND cp maffilter.1.txt maffilter.1
      COMMAND ${DOC_COMPRESS} -f maffilter.1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
      )
    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "man/maffilter.1.${DOC_COMPRESS_EXT}")
  ELSE()
    MESSAGE(FATAL_ERROR "Program nroff required but not found.")
  ENDIF()
ENDIF(MAN)

IF(INFO)
  FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
  IF(MAKEINFO_EXE)
    MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
    MESSAGE("   Adding targets: info, html")
 
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
    
    ADD_CUSTOM_TARGET(info
      ${ADD_INFO_TO}
      COMMAND ${MAKEINFO_EXE} maffilter.texi
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
      )

    ADD_CUSTOM_TARGET(html
      COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely maffilter.texi
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
      )
  ELSE(MAKEINFO_EXE)
    MESSAGE(FATAL_ERROR p"Program makeinfo required but not found.")
  ENDIF(MAKEINFO_EXE)
ENDIF(INFO)

FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)
IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} maffilter.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "doc/maffilter.pdf;doc/maffilter.aux;doc/maffilter.cp;doc/maffilter.fn;doc/maffilter.info;doc/maffilter.ky;doc/maffilter.log;doc/maffilter.pg;doc/maffilter.toc;doc/maffilter.tp;doc/maffilter.vr")
ENDIF(PDFTEX_EXE)

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  IF(BUILD_STATIC)
    SET(${lib_name}_STATIC_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_STATIC_LIBRARY NAMES ${${lib_name}_STATIC_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_STATIC_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
  ELSE()
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  ENDIF()
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

#Find the libraries. The order is very important for static linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl-omics Bpp/Seq/Io/Maf/AbstractDistanceEstimationMafIterator.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq-omics Bpp/Seq/Feature/SequenceFeature.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

#Find boost libraries
SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.34.0 COMPONENTS iostreams )

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost found: ${Boost_INCLUDE_DIRS}")
ENDIF()

# Find the zlib installation
IMPROVED_FIND_LIBRARY(LIBS z zlib.h)
#This does not work with static linkage!!!
#FIND_PACKAGE(ZLIB REQUIRED)
#IF(ZLIB_FOUND)
#  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
#  SET(LIBS ${LIBS} ${ZLIB_LIBRARIES})
#  MESSAGE("-- Zlib found here:")
#  MESSAGE("   includes: ${ZLIB_INCLUDE_DIR}")
#  MESSAGE("   libraries: ${ZLIB_LIBRARIES}")
#ENDIF()

# Find the bz2 installation
IMPROVED_FIND_LIBRARY(LIBS bz2 bzlib.h)
#This does not work with static linkage!!!
#FIND_PACKAGE(BZip2 REQUIRED)
#IF(BZIP2_FOUND)
#  INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
#  SET(LIBS ${LIBS} ${BZIP2_LIBRARIES})
#  MESSAGE("-- BZip2 found here:")
#  MESSAGE("   includes: ${BZIP2_INCLUDE_DIR}")
#  MESSAGE("   libraries: ${BZIP2_LIBRARIES}")
#ENDIF()

# Subdirectories
ADD_SUBDIRECTORY(MafFilter)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(man)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "maffilter")
SET(CPACK_PACKAGE_VENDOR "Julien Dutheil")
SET(CPACK_PACKAGE_VERSION "1.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Filtering of genome alignment in the Multiple Alignment Format (MAF)")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.git"
 ".*\\\\.gz"
 ".*\\\\.zip"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 ".*\\\\.dmg"
 ".*\\\\.sh"
 ".*\\\\..*\\\\.swp"
 ".*stamp"
 "\\\\.sh"
 "\\\\..*\\\\.swp"
 "MafFilter/\\\\..*"
 "man/.*\\\\.1.${DOC_COMPRESS_EXT}"
 "debian/tmp"
 "debian/maffilter/"
 "debian/maffilter\\\\.substvars"
 "debian/maffilter\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "install_manifest.txt"
 "DartConfiguration.tcl"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
# 'clean' is not (yet) a first class target. However, we need to clean the directories before building the sources:
IF("${CMAKE_GENERATOR}" MATCHES "Make")
  ADD_CUSTOM_TARGET(make_clean
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )
  ADD_DEPENDENCIES(dist make_clean)

  ADD_CUSTOM_TARGET(make_clean_man
    COMMAND rm -f *.${DOC_COMPRESS_EXT}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
  )
  ADD_DEPENDENCIES(dist make_clean_man)
ENDIF()

IF(NOT NO_DEP_CHECK)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(origdist COMMAND cp ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz ../${CPACK_DEBSOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(origdist dist)
ADD_CUSTOM_TARGET(deb dpkg-buildpackage -uc -us -i${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(deb origdist)
ADD_DEPENDENCIES(deb info man)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

ENDIF(UNIX)

ENDIF(NOT NO_DEP_CHECK)
