# SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(
    maffilter
    MafFilter.cpp
    OutputAsFeaturesMafIterator.cpp
    SystemCallMafIterator.cpp
    TreeBuildingSystemCallMafIterator.cpp
    PopGenMafStatistics.cpp
)
set(maffilter-targets maffilter)

foreach(target ${maffilter-targets})
    # Link (static or shared)
    if(BUILD_STATIC)
	    target_link_libraries(${target} ${BPP_LIBS_STATIC} ${LIBS})
        set_target_properties(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
	target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
    else(BUILD_STATIC)
	    target_link_libraries(${target} ${BPP_LIBS_SHARED} ${LIBS})
        set_target_properties(
            ${target}
            PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        )
    endif(BUILD_STATIC)
endforeach(target)

install(TARGETS ${maffilter-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
