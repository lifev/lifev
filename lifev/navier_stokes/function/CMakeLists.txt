SET(function_HEADERS
function/KimMoin.hpp
function/RossEthierSteinmanDec.hpp
function/RossEthierSteinmanInc.hpp
function/Womersley.hpp
CACHE INTERNAL "")

SET(function_SOURCES
function/KimMoin.cpp
function/RossEthierSteinmanDec.cpp
function/RossEthierSteinmanInc.cpp
function/Womersley.cpp
CACHE INTERNAL "")

ADD_SUBDIRECTORY(bessel)
