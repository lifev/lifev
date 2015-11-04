SET( solver_IonicModels_HEADERS
  solver/IonicModels/IonicAlievPanfilov.hpp
  solver/IonicModels/IonicFitzHughNagumo.hpp
  solver/IonicModels/IonicLuoRudyI.hpp
  solver/IonicModels/IonicMitchellSchaeffer.hpp
  solver/IonicModels/IonicTenTusscher06.hpp
  solver/IonicModels/IonicMinimalModel.hpp
  solver/IonicModels/ElectroIonicModel.hpp
  solver/IonicModels/IonicFox.hpp
  solver/IonicModels/IonicHodgkinHuxley.hpp
  solver/IonicModels/IonicNoblePurkinje.hpp
  solver/IonicModels/IonicGoldbeter.hpp
CACHE INTERNAL "")

SET( solver_IonicModels_SOURCES
  solver/IonicModels/IonicAlievPanfilov.cpp
  solver/IonicModels/IonicFitzHughNagumo.cpp
  solver/IonicModels/IonicLuoRudyI.cpp
  solver/IonicModels/IonicMitchellSchaeffer.cpp
  solver/IonicModels/IonicTenTusscher06.cpp
  solver/IonicModels/IonicMinimalModel.cpp
  solver/IonicModels/IonicFox.cpp
  solver/IonicModels/ElectroIonicModel.cpp
  solver/IonicModels/IonicHodgkinHuxley.cpp
  solver/IonicModels/IonicNoblePurkinje.cpp
  solver/IonicModels/IonicGoldbeter.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${solver_IonicModels_HEADERS})
  STRING(REPLACE "solver/IonicModels/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/lifev/electrophysiology/solver/IonicModels"
  COMPONENT ${PACKAGE_NAME}
)
