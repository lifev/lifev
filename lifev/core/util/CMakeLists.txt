SET(util_HEADERS
util/CBlas.hpp
util/CLapack.hpp
util/Displayer.hpp
util/EncoderBase64.hpp
util/Factory.hpp
util/FactoryPolicy.hpp
util/FactorySingleton.hpp
util/FactoryTypeInfo.hpp
util/FortranWrapper.hpp
util/LifeAssert.hpp
util/LifeAssertSmart.hpp
util/LifeChrono.hpp
util/LifeDebug.hpp
util/ParserDefinitions.hpp
util/Parser.hpp
util/ParserSpiritGrammar.hpp
util/StringData.hpp
util/StringUtility.hpp
util/Switch.hpp
CACHE INTERNAL "")

SET(util_SOURCES
util/Displayer.cpp
util/EncoderBase64.cpp
util/FactoryTypeInfo.cpp
util/LifeAssertSmart.cpp
util/LifeDebug.cpp
util/Parser.cpp
util/StringData.cpp
util/StringUtility.cpp
util/Switch.cpp
CACHE INTERNAL "")
