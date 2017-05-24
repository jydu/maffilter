# CMake script for MafFilter package
# Author: Julien Dutheil
# Created: 27/04/2010

ADD_EXECUTABLE(maffilter MafFilter.cpp OutputAsFeaturesMafIterator.cpp SystemCallMafIterator.cpp TreeBuildingSystemCallMafIterator.cpp)
set (maffilter-targets
  maffilter
  )

foreach (target ${maffilter-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    target_link_libraries (${target} ${LIBS})
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
    target_link_libraries (${target} ${LIBS})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${maffilter-targets} DESTINATION bin)
