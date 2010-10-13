SET(FEMRI_SYMBOLIC_SRCS
  vtkfemriUnstructuredGridFourierIntegrator.cxx
  vtkfemriUnstructuredGridKSpaceGenerator.cxx
  )

#IF(BUILD_TESTING)
#  ENABLE_TESTING()
#  ADD_EXECUTABLE(femrisymbolictesting femriSymbolicTesting.cxx ${FEMRI_SYMBOLIC_SRCS})
#  TARGET_LINK_LIBRARIES(femrisymbolictesting ginac cln stdc++ vtkCommon vtkIO)
#  ADD_TEST(femrisymbolictest femrisymbolictesting)
#ENDIF(BUILD_TESTING)

CONFIGURE_FILE(
  ${FEMRI_SOURCE_DIR}/femriSymbolic/vtkfemriSymbolicConfigure.h.in
  ${FEMRI_BINARY_DIR}/femriSymbolic/vtkfemriSymbolicConfigure.h
  )

FIND_PATH(GINAC_INCLUDE_DIR ginac.h /usr/include/ginac /usr/local/include/ginac)
FIND_PATH(GINAC_LIB_DIR libginac.so /usr/lib /usr/local/lib)
INCLUDE_DIRECTORIES(${GINAC_INCLUDE_DIR})
LINK_DIRECTORIES(${GINAC_LIB_DIR})

ADD_LIBRARY(femriSymbolic ${FEMRI_SYMBOLIC_SRCS})
TARGET_LINK_LIBRARIES(femriSymbolic ginac cln stdc++ vtkCommon vtkIO vtkImaging)

#ADD_EXECUTABLE(femrisymbolicapp femriSymbolicApp.cxx ${FEMRI_SYMBOLIC_SRCS})
#TARGET_LINK_LIBRARIES(femrisymbolicapp ginac cln stdc++ vtkCommon vtkIO vtkImaging)

INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriSymbolic)
INSTALL_FILES(${FEMRI_INSTALL_INCLUDE_DIR} .h ${FEMRI_SYMBOLIC_SRCS})

IF (FEMRI_WRAP_PYTHON)
  CONFIGURE_FILE(${FEMRI_SOURCE_DIR}/femriSymbolic/femrisymbolic.py ${FEMRI_BINARY_DIR}/femriSymbolic/femrisymbolic.py COPYONLY)
  INSTALL_FILES(${FEMRI_INSTALL_LIB_DIR} FILES ${FEMRI_BINARY_DIR}/femriSymbolic/femrisymbolic.py)
  VTK_WRAP_PYTHON3(femriSymbolicPython FEMRI_SYMBOLIC_PYTHON_SRCS "${FEMRI_SYMBOLIC_SRCS}")
  ADD_LIBRARY(femriSymbolicPythonD ${FEMRI_SYMBOLIC_PYTHON_SRCS})
  ADD_LIBRARY(femriSymbolicPython MODULE femriSymbolicPythonInit.cxx)
  TARGET_LINK_LIBRARIES(femriSymbolicPythonD femriSymbolic vtkCommon vtkCommonPythonD vtkIO vtkIOPythonD vtkImaging vtkImagingPythonD)
  TARGET_LINK_LIBRARIES (femriSymbolicPython femriSymbolicPythonD)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriSymbolicPythonD)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriSymbolicPython)
ENDIF (FEMRI_WRAP_PYTHON)

