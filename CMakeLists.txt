# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Global parameters
cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(maffilter CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile options
add_compile_options(-Wall -Weffc++ -Wshadow -Wconversion)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif()


#static linkage?
if(NOT BUILD_STATIC)
    set(BUILD_STATIC FALSE CACHE BOOL "Enable static linkage." FORCE)
endif()
if(BUILD_STATIC)
    message(STATUS "Static linkage requested.")
endif()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if(NOT COMPRESS_PROGRAM)
    set(COMPRESS_PROGRAM
        gzip
        CACHE STRING
        "Set program for compressing documentation"
        FORCE
    )
endif()
find_program(
    COMPRESS_BIN
    NAMES ${COMPRESS_PROGRAM}
    DOC "${COMPRESS_PROGRAM} compression program"
)
if(NOT COMPRESS_BIN)
    message(
        STATUS
        "${COMPRESS_PROGRAM} program not found, text doc will not be compressed"
    )
else()
    # Deduce COMPRESS_EXT for known compression programs if not set
    if(NOT COMPRESS_EXT)
        if(${COMPRESS_PROGRAM} STREQUAL "gzip")
            set(COMPRESS_EXT "gz")
        elseif(${COMPRESS_PROGRAM} STREQUAL "bzip2")
            set(COMPRESS_EXT "bz2")
        else()
            set(COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
        endif()
    endif()
    # Generate command line args (always add -c to output compressed file to stdout)
    if(${COMPRESS_PROGRAM} STREQUAL "gzip")
        # -n for no timestamp in files (reproducible builds)
        # -9 for maximum compression (lintian error)
        set(COMPRESS_ARGS -c -n -9)
    else()
        set(COMPRESS_ARGS -c)
    endif()
    message(
        STATUS
        "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}"
    )
endif()

# Find dependencies (add install directory to search)
if(CMAKE_INSTALL_PREFIX)
    set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif(CMAKE_INSTALL_PREFIX)

include(GNUInstallDirs)
find_package(bpp-phyl-omics 3.0.0 REQUIRED)
find_package(bpp-popgen 3.0.0 REQUIRED)

#Find boost libraries
set(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
set(Boost_USE_MULTITHREADED ON)
find_package(Boost CONFIG REQUIRED)
#find_package(Boost 1.34.0 COMPONENTS iostreams)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${Boost_LIBRARIES})
    message("-- Boost libraries found here:")
    message("   includes: ${Boost_INCLUDE_DIRS}")
    message("   libraries: ${Boost_LIBRARIES}")
endif()

find_package(Eigen3 3.3 REQUIRED PATHS /usr/lib NO_MODULE)

#here is a useful function:
macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
    #start:
    find_path(${lib_name}_INCLUDE_DIR ${include_to_find})
    include_directories(${${lib_name}_INCLUDE_DIR})

    if(BUILD_STATIC)
        set(${lib_name}_STATIC_NAMES lib${lib_name}.a)
        find_library(
            ${lib_name}_STATIC_LIBRARY
            NAMES ${${lib_name}_STATIC_NAMES}
            PATH_SUFFIXES lib${LIB_SUFFIX}
        )
        if(${lib_name}_STATIC_LIBRARY)
            message("-- Library ${lib_name} found here:")
            message("   includes: ${${lib_name}_INCLUDE_DIR}")
            message("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
        else()
            message(FATAL_ERROR "${lib_name} required but not found.")
        endif()
        #add the dependency:
        set(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
    else()
        set(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
        find_library(
            ${lib_name}_LIBRARY
            NAMES ${${lib_name}_NAMES}
            PATH_SUFFIXES lib${LIB_SUFFIX}
        )
        if(${lib_name}_LIBRARY)
            message("-- Library ${lib_name} found here:")
            message("   includes: ${${lib_name}_INCLUDE_DIR}")
            message("   dynamic libraries: ${${lib_name}_LIBRARY}")
        else()
            message(FATAL_ERROR "${lib_name} required but not found.")
        endif()
        #add the dependency:
        set(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
    endif()
endmacro(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
if(CMAKE_INSTALL_PREFIX)
    set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif(CMAKE_INSTALL_PREFIX)

# Find the zlib installation
improved_find_library(LIBS z zlib.h)
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
improved_find_library(LIBS bz2 bzlib.h)
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
add_subdirectory(MafFilter)
add_subdirectory(doc)
add_subdirectory(man)

# Packager
set(CPACK_PACKAGE_NAME "maffilter")
set(CPACK_PACKAGE_VENDOR "Julien Y. Dutheil")
set(CPACK_PACKAGE_VERSION "1.4.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "4")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Filtering of genome alignment in the Multiple Alignment Format (MAF)"
)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
set(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL")
set(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/examples/"
    "/\\\\.git/"
    "/\\\\.gitignore"
    ${CPACK_SOURCE_IGNORE_FILES}
)

set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)
set(CPACK_DEBSOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig"
)
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

if(UNIX)
    #This creates rpm packages:
    add_custom_target(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
    add_dependencies(rpm dist)
endif(UNIX)
