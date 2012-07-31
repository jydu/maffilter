# CMake script for GenomeTools package
# Author: Julien Dutheil
# Created: 27/04/2010

ADD_EXECUTABLE(mafinfo MafInfo.cpp)
TARGET_LINK_LIBRARIES(mafinfo ${LIBS} ${Boost_LIBRARIES})

ADD_EXECUTABLE(mafqdist MafQDist.cpp)
TARGET_LINK_LIBRARIES(mafqdist ${LIBS} ${Boost_LIBRARIES})

ADD_EXECUTABLE(maffilter MafFilter.cpp MafIterators.cpp)
TARGET_LINK_LIBRARIES(maffilter ${LIBS} ${Boost_LIBRARIES})

# Install progs
INSTALL(TARGETS
  mafinfo
  mafqdist
  maffilter
  DESTINATION bin)

