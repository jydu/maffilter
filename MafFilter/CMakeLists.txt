# CMake script for MafFilter package
# Author: Julien Dutheil
# Created: 27/04/2010

ADD_EXECUTABLE(maffilter MafFilter.cpp OutputAsFeaturesMafIterator.cpp)
TARGET_LINK_LIBRARIES(maffilter ${LIBS} ${Boost_LIBRARIES})
SET_TARGET_PROPERTIES(maffilter PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})

# Install progs
INSTALL(TARGETS
  maffilter
  DESTINATION bin)

