SET(util_HEADERS
  util/FactoryTypeInfo.hpp
  util/Switch.hpp
  util/EncoderBase64.hpp
  util/LifeDebug.hpp
  util/StringUtility.hpp
  util/LifeAssertSmart.hpp
  util/Parser.hpp
  util/FactorySingleton.hpp
  util/StringData.hpp
  util/ParserSpiritGrammar.hpp
  util/FortranWrapper.hpp
  util/Factory.hpp
  util/LifeAssert.hpp
  util/LifeChrono.hpp
  util/Displayer.hpp
  util/ParserDefinitions.hpp
  util/FactoryPolicy.hpp
  util/WallClock.hpp
CACHE INTERNAL "")

SET(util_SOURCES
  util/LifeDebug.cpp
  util/StringData.cpp
  util/EncoderBase64.cpp
  util/StringUtility.cpp
  util/LifeAssertSmart.cpp
  util/Switch.cpp
  util/Parser.cpp
  util/FactoryTypeInfo.cpp
  util/Displayer.cpp
  util/WallClock.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${util_HEADERS})
  STRING(REPLACE "util/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/core/util"
  COMPONENT ${PACKAGE_NAME}
)
