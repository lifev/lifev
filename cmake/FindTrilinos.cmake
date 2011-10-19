# This three lines are left to be similar to Radu, but I do not understand what they do.
set(Trilinos_FOUND TRUE)
set(Trilinos_ERROR_REASONS)
set(Trilinos_ALREADY_CACHED FALSE)

# Here I am looking for TrilinosConfig.cmake and I will import it. 
FIND_PACKAGE(Trilinos PATHS ${Trilinos_ROOT}/include)

# Stop cmake if Trilinos is not found.
IF(NOT Trilinos_FOUND)
  MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()

# Here it will be better just to raise a warning or have if(USE_TRILINOS_COMPILERS)
# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

# Set all the include directories to compile against Trilinos and its dependencies
SET (Trilinos_INCLUDE_DIR  ${Trilinos_DIR} ${Trilinos_TPL_INCLUDE_DIRS})

# Now we created the list of all the Trilinos libraries and TPL we need to link agaist to
#First I define the directory where Trilinos libs are
SET( tmp  "-L${Trilinos_ROOT}/lib")	
# Then I concatenate all the Trilinos libraries -lXXXX
FOREACH(LIB_FILE ${Trilinos_LIBRARIES})
	SET (tmp "${tmp} -l${LIB_FILE}")
ENDFOREACH(LIB_FILE ${Trilinos_LIBRARIES})    
#Finally I set Trilinos_LIBRARIES to be the set of Trilinos libraries plus the other TPL of Trilinos
SET (Trilinos_LIBRARIES ${tmp} ${Trilinos_TPL_LIBRARIES})

#Finally I check that all the Trilinos package we need are present.
#I define all the needed packages
SET( REQUIRED_PACKAGES "amesos;anasazi;aztecoo;belos;epetra;epetraext;ifpack;ml;teuchos;zoltan")
#Now I check if all the required packages are in the list of Trilinos Packages. (This might be done better with a find...)
FOREACH( PACKAGE ${REQUIRED_PACKAGES})
	STRING(TOUPPER ${PACKAGE} UPACK)
	FOREACH(Tpack ${Trilinos_PACKAGE_LIST})
		STRING(TOUPPER ${Tpack} TPACK)
		IF(${UPACK} STREQUAL ${TPACK})
			SET(${UPACK}_FOUND TRUE)
		ENDIF()
	ENDFOREACH(Tpack ${Trilinos_PACKAGE_LIST})
	IF(${UPACK}_FOUND)
		set( HAVE_TRILINOS_${UPACK} TRUE)
	ELSE()
		MESSAGE( FATAL_ERROR "Could not find ${PACKAGE}")
	ENDIF()
ENDFOREACH(PACKAGE ${REQUIRED_PACKAGES})

#If XXX is already linked through Trilinos. If so I set the variable XXX_IS_IN_TRILINOS to true
FOREACH(tpl ${Trilinos_TPL_LIST})
	STRING(TOUPPER ${tpl} TPL)
	SET(${TPL}_IS_IN_TRILINOS TRUE)
	MESSAGE( "Library ${TPL} already included through Trilinos" )
ENDFOREACH(tpl ${Trilinos_TPL_LIST})
