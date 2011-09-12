add_library (lifecore Displayer.cpp EncoderBase64.cpp FactoryTypeInfo.cpp LifeAssertSmart.cpp LifeDebug.cpp LifeV.cpp LifeVersion.cpp Parser.cpp StringData.cpp StringUtility.cpp Switch.cpp)
install (TARGETS lifecore DESTINATION lib)
install (FILES CBlas.hpp CLapack.hpp Displayer.hpp EncoderBase64.hpp Factory.hpp FactoryPolicy.hpp FactorySingleton.hpp FactoryTypeInfo.hpp FortranWrapper.hpp LifeAssert.hpp LifeAssertSmart.hpp LifeChrono.hpp LifeDebug.hpp LifeVersion.hpp LifeV.hpp ParserDefinitions.hpp Parser.hpp ParserSpiritGrammar.hpp StringData.hpp StringUtility.hpp Switch.hpp DESTINATION include/life/lifecore)

