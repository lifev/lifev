###################################################################################################
#
#                       This file is part of the LifeV Applications                        
#                Copyright (C) 2001-(>>>YEAR<<<) EPFL, Politecnico di Milano, INRIA          
#
#      Author(s): (>>>USER_NAME<<<) <(>>>AUTHOR<<<)>
#           Date: (>>>DATE<<<)
#  License Terms: GNU GPL
#
###################################################################################################

AT_SETUP([test_NameOfTheTest])    #Name of the test
AT_KEYWORDS([])
AT_DATA([data.txt],               #Data file
[[

(>>>POINT<<<)
# Place here the content of the GetPot data file for the night test.

]])
AT_CHECK([ln -sf ../../data/mesh/inria/Mesh &&
		  mpirun -n 1 ../../test_TemplateTest/test_TemplateTest -c],[0],[ignore],[ignore])
AT_CLEANUP([FilesCreatedByTheTest1.txt 
            FilesCreatedByTheTest2.txt 
            ...
          ]) # Files to be removed at the end of the test
