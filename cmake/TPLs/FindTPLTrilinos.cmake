INCLUDE(TribitsTplDeclareLibraries)

# TRIBITS_TPL_DECLARE_LIBRARIES( Trilinos
#   REQUIRED_HEADERS Epetra_Comm.h
#   REQUIRED_LIBS_NAMES "epetra"
#   )

if(Trilinos_LIBRARY_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} "${Trilinos_LIBRARY_DIRS}/cmake/Trilinos")
endif()
if (Trilinos_INCLUDE_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} ${Trilinos_INCLUDE_DIRS})
endif()

MESSAGE("Trilinos_DIR: ${Trilinos_DIR}")

# Here I am looking for TrilinosConfig.cmake and I will import it. 
find_package (Trilinos NO_MODULE HINTS ${Trilinos_DIR})

# Stop cmake if Trilinos is not found.
if (NOT Trilinos_FOUND)
  message (FATAL_ERROR "Could not find Trilinos!")
endif ()

if("${Trilinos_VERSION_MINOR}" GREATER 6)
  set (HAVE_TRILINOS_GT_10_6 TRUE)
  message (STATUS "Using Trilinos > 10.6 " ${Trilinos_VERSION_MINOR})
else()
  message (STATUS "Using Trilinos <= 10.6 " ${Trilinos_VERSION_MINOR})
endif ()

# Here it will be better just to raise a warning or have if(USE_TRILINOS_COMPILERS)
# Make sure to use same compilers and flags as Trilinos
IF(NOT ${CMAKE_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos CXX compiler")
  MESSAGE(STATUS "CMAKE_CXX_COMPILER:    " ${CMAKE_CXX_COMPILER})
  MESSAGE(STATUS "Trilinos_CXX_COMPILER: " ${Trilinos_CXX_COMPILER})
ENDIF()
IF(NOT ${CMAKE_C_COMPILER} STREQUAL ${Trilinos_C_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos C compiler")
  MESSAGE(STATUS "CMAKE_C_COMPILER:    " ${CMAKE_C_COMPILER})
  MESSAGE(STATUS "Trilinos_C_COMPILER: " ${Trilinos_C_COMPILER})
ENDIF()
IF(NOT ${CMAKE_Fortran_COMPILER} STREQUAL ${Trilinos_Fortran_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos Fortran compiler")
  MESSAGE(STATUS "CMAKE_Fortran_COMPILER:    " ${CMAKE_C_COMPILER})
  MESSAGE(STATUS "Trilinos_Fortran_COMPILER: " ${Trilinos_C_COMPILER})
ENDIF()

# Optional Packages (to be moved outside with COMPONENTS ...)
list (APPEND LifeV_OPTIONAL_Trilinos_PKGS
  "NOX" "Thyra" "Rythmos" "Teko" "Stratimikos" "Isorropia" "ShyLU")

# Required packages (to be moved outside, like REQUIRED COMPONENTS ...)
list (APPEND LifeV_REQUIRED_Trilinos_PKGS
  "ML" "Ifpack" "Amesos" "Anasazi" "Belos" "AztecOO" "Zoltan" "EpetraExt" "Epetra" "Teuchos")

# Start scanning Trilinos configuration
foreach (TYPE IN ITEMS "OPTIONAL" "REQUIRED")
  foreach (PKG IN LISTS LifeV_${TYPE}_Trilinos_PKGS)
    # Look for PKG
    list (FIND Trilinos_PACKAGE_LIST "${PKG}" PKG_FOUND)
    if (PKG_FOUND GREATER -1)
      # Found! Let's announce it!
      message (STATUS "Trilinos :: ${PKG} Found!")
      list (APPEND LifeV_Trilinos_LIBRARIES "${${PKG}_LIBRARIES}")
      list (APPEND LifeV_Trilinos_TPL_INCLUDE_DIRS "${${PKG}_TPL_INCLUDE_DIRS}")
      list (APPEND LifeV_Trilinos_TPL_LIST "${${PKG}_TPL_LIST}")
      string (TOUPPER ${PKG} UPKG)
      set (${UPKG}_FOUND True)
      set (HAVE_TRILINOS_${UPKG} True)
    else ()
      if (TYPE STREQUAL "REQUIRED")
        message (FATAL_ERROR "Trilinos :: ${PKG} NOT Found!")
      else ()
        message (WARNING "Trilinos :: ${PKG} NOT Found! Some test might not compile properly ...")
      endif ()
    endif ()
  endforeach (PKG)
endforeach (TYPE)

# Cleaning duplicates
list (REVERSE LifeV_Trilinos_TPL_LIST)
list (REMOVE_DUPLICATES LifeV_Trilinos_TPL_LIST)
list (REVERSE LifeV_Trilinos_TPL_LIST)
list (REVERSE LifeV_Trilinos_LIBRARIES)
list (REMOVE_DUPLICATES LifeV_Trilinos_LIBRARIES)
list (REVERSE LifeV_Trilinos_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
list (REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
set (LifeV_Trilinos_TPL_LIBRARIES ${Trilinos_TPL_LIBRARIES})

list (REMOVE_DUPLICATES LifeV_Trilinos_TPL_INCLUDE_DIRS)

list (APPEND LifeV_Trilinos_INCLUDE_DIRS
  ${Trilinos_INCLUDE_DIRS}
  ${LifeV_Trilinos_TPL_INCLUDE_DIRS})
# I think there's a better way to handle this ... CMake
# should take care of -L or -l or -rpath ...
set (LifeV_Trilinos_LIBS "-L${Trilinos_LIBRARY_DIRS}")
foreach (LIB IN LISTS LifeV_Trilinos_LIBRARIES)
  set (LifeV_Trilinos_LIBS "${LifeV_Trilinos_LIBS} -l${LIB}")
endforeach (LIB)
set (LifeV_Trilinos_LIBS ${LifeV_Trilinos_LIBS} ${LifeV_Trilinos_TPL_LIBRARIES})

# TPLs
foreach (TPL IN ITEMS "ParMETIS" "Boost" "LAPACK" "BLAS" "UMFPACK" "SuperLU" "SuperLUDist" "HDF5")
    list (FIND LifeV_Trilinos_TPL_LIST ${TPL} TPL_FOUND)
  if (TPL_FOUND GREATER -1)
    string (TOUPPER ${TPL} UTPL)
    set (${UTPL}_IS_IN_TRILINOS True)
  endif()
endforeach (TPL)

# Filling variables needed by the TriBITS system
set (TPL_Trilinos_INCLUDE_DIRS ${LifeV_Trilinos_INCLUDE_DIRS})
set (TPL_Trilinos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS})
set (TPL_Trilinos_LIBRARIES ${LifeV_Trilinos_LIBS})
