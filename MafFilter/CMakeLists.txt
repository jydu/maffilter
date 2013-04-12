# CMake script for MafFilter package
# Author: Julien Dutheil
# Created: 27/04/2010

ADD_EXECUTABLE(mafqdist MafQDist.cpp)
TARGET_LINK_LIBRARIES(mafqdist ${LIBS} ${Boost_LIBRARIES})

ADD_EXECUTABLE(maffilter MafFilter.cpp OutputAsFeaturesMafIterator.cpp)
TARGET_LINK_LIBRARIES(maffilter ${LIBS} ${Boost_LIBRARIES})

# Install progs
INSTALL(TARGETS
  mafqdist
  maffilter
  DESTINATION bin)

