SET(stimulus_HEADERS
  stimulus/ElectroStimulus.hpp  
  stimulus/StimulusPMJ.hpp
  stimulus/StimulusSingleSource.hpp
  stimulus/StimulusPacingProtocol.hpp
CACHE INTERNAL "")

SET(stimulus_SOURCES
  stimulus/ElectroStimulus.cpp  
  stimulus/StimulusPMJ.cpp
  stimulus/StimulusSingleSource.cpp
  stimulus/StimulusPacingProtocol.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${stimulus_HEADERS})
  STRING(REPLACE "stimulus/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION
"${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/stimulus"
  COMPONENT ${PACKAGE_NAME}
)
