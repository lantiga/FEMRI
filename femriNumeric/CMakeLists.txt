
SET(FEMRI_NUMERIC_SRCS
#  vtkfemriUnstructuredGridMRSlicer.cxx
#  vtkSubdivideSelectTetra.cxx
#  vtkSubdivideSelectQuadraticTetra.cxx
  vtkfemriUnstructuredGridNumericKSpaceGenerator.cxx
  vtkfemriPolyDataNumericKSpaceGenerator.cxx
  vtkfemriGaussQuadrature.cxx
  vtkfemriOptimalQuadratureOrderCalculator.cxx
  femriIntegrationDesignCurve.cxx
  femriIntegrationDesignCurveInterpolator.cxx
  )

SET(FEMRI_NUMERIC_TEST_SRCS
  femriInterpolatorTest.cxx
  )

SET_SOURCE_FILES_PROPERTIES(
  femriIntegrationDesignCurve.cxx
  femriIntegrationDesignCurveInterpolator.cxx
  WRAP_EXCLUDE)

CONFIGURE_FILE(
  ${FEMRI_SOURCE_DIR}/femriNumeric/vtkfemriNumericConfigure.h.in
  ${FEMRI_BINARY_DIR}/femriNumeric/vtkfemriNumericConfigure.h
  )

ADD_LIBRARY(femriNumeric ${FEMRI_NUMERIC_SRCS})
TARGET_LINK_LIBRARIES(femriNumeric femriCommon vtkCommon vtkFiltering vtkIO vtkImaging vtkGraphics)

ADD_EXECUTABLE(femriNumericInterpolatorTest ${FEMRI_NUMERIC_TEST_SRCS})
TARGET_LINK_LIBRARIES(femriNumericInterpolatorTest femriNumeric)

INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriNumeric)
INSTALL_FILES(${FEMRI_INSTALL_INCLUDE_DIR} .h ${FEMRI_NUMERIC_SRCS})

IF (FEMRI_WRAP_PYTHON)
  CONFIGURE_FILE(${FEMRI_SOURCE_DIR}/femriNumeric/femrinumeric.py ${FEMRI_BINARY_DIR}/femriNumeric/femrinumeric.py COPYONLY)
  INSTALL_FILES(${FEMRI_INSTALL_LIB_DIR} FILES ${FEMRI_BINARY_DIR}/femriNumeric/femrinumeric.py)
  VTK_WRAP_PYTHON3(femriNumericPython FEMRI_NUMERIC_PYTHON_SRCS "${FEMRI_NUMERIC_SRCS}")
  ADD_LIBRARY(femriNumericPythonD ${FEMRI_NUMERIC_PYTHON_SRCS})
  ADD_LIBRARY(femriNumericPython MODULE femriNumericPythonInit.cxx)
  TARGET_LINK_LIBRARIES(femriNumericPythonD femriNumeric femriCommon femriCommonPythonD vtkCommon vtkCommonPythonD vtkIO vtkIOPythonD vtkImaging vtkImagingPythonD vtkGraphics vtkGraphicsPythonD)
  TARGET_LINK_LIBRARIES (femriNumericPython femriNumericPythonD)
  IF(WIN32 AND NOT CYGWIN)
     SET_TARGET_PROPERTIES(femriNumericPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriNumericPythonD)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriNumericPython)
ENDIF (FEMRI_WRAP_PYTHON)

