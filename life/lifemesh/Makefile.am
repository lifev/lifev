
SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc .cc

INCLUDES 	= -I$(top_srcdir)/life/lifecore -I$(top_srcdir)/life/lifefem
AM_CPPFLAGS	= -DTHREEDIM


include_HEADERS		= bareItems.hpp basisElSh.hpp bcCond.hpp bcVector.hpp \
			  dof.hpp dofInterface3Dto2D.hpp dofInterface3Dto3D.hpp dofInterfaceBase.hpp \
			  dofInterfaceHandler.hpp geo0D.hpp geoElement.hpp geoND.hpp \
                          identifier.hpp markers.hpp markers_base.hpp meshEntity.hpp \
                          meshMotion.hpp mesh_util.hpp mesh_util_base.hpp readMesh2D.hpp \
                          readMesh3D.hpp regionMesh2D.hpp regionMesh3D.hpp 

lib_LTLIBRARIES		= liblifemesh.la
liblifemesh_la_SOURCES	= bareItems.cpp basisElSh.cpp bcCond.cpp bcVector.cpp \
                          dof.cpp dofInterface3Dto2D.cpp dofInterface3Dto3D.cpp dofInterfaceBase.cpp \
                          dofInterfaceHandler.cpp geo0D.cpp identifier.cpp markers_base.cpp \
                          meshMotion.cpp mesh_util.cpp mesh_util_base.cpp readMesh3D.cpp \
                          regionMesh3D.cpp 

liblifemesh_la_LIBADD	= ../lifecore/liblifecore.la
