TESTSUITE README:

This folder contains the template files for creating a new test in the testsuite.

Suppose that we want to create a new test with the name "MyTest", than we have to:
1) create a new folder with the name "test_MyTest" in the testsuite main directory;
2) copy/paste all the files from "test_TemplateTest" to "test_MyTest";
3) modify the main.cpp and Makefile.am for running the new test;
4) modify the data.txt and testsuite.at, adding the parameters;
5) check that the test is correctly executed when typing "make check".
   To do it edit the following files:
   
   testsuite/Makefile.am
   testsuite/testsuite.at
   
6) add the test to the global doxygen documentation by inserting 

   $(top_srcdir)/testsuite/test_MyTest

   in the INPUT section of Doxyfile.am  

For any additional information about this procedure, please read the "LifeV Development Guidelines"
or contact the LifeV administrator.