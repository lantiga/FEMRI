SET(FEMRI_COMMON_SRCS
  vtkfemriKSpaceZeroPadder.cxx
  vtkfemriKSpaceGenerator.cxx
  vtkfemriKSpaceShift.cxx
  vtkfemri2DBoxKSpaceGenerator.cxx
  vtkfemri2DCylinderKSpaceGenerator.cxx
  vtkfemri3DBoxKSpaceGenerator.cxx
  vtkfemri3DCylinderKSpaceGenerator.cxx
  vtkfemri3DSphereKSpaceGenerator.cxx
  vtkfemriBoundaryIntegralKSpaceGenerator.cxx
  )

CONFIGURE_FILE(
  ${FEMRI_SOURCE_DIR}/femriCommon/vtkfemriCommonConfigure.h.in
  ${FEMRI_BINARY_DIR}/femriCommon/vtkfemriCommonConfigure.h
  ) 

ADD_LIBRARY(femriCommon ${FEMRI_COMMON_SRCS})
TARGET_LINK_LIBRARIES(femriCommon vtkCommon vtkImaging)

INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriCommon)
INSTALL_FILES(${FEMRI_INSTALL_INCLUDE_DIR} .h ${FEMRI_COMMON_SRCS})

IF (FEMRI_WRAP_PYTHON)
  CONFIGURE_FILE(${FEMRI_SOURCE_DIR}/femriCommon/femricommon.py ${FEMRI_BINARY_DIR}/femriCommon/femricommon.py COPYONLY)
  INSTALL_FILES(${FEMRI_INSTALL_LIB_DIR} FILES ${FEMRI_BINARY_DIR}/femriCommon/femricommon.py)
  VTK_WRAP_PYTHON3(femriCommonPython FEMRI_COMMON_PYTHON_SRCS "${FEMRI_COMMON_SRCS}")
  ADD_LIBRARY(femriCommonPythonD ${FEMRI_COMMON_PYTHON_SRCS})
  ADD_LIBRARY(femriCommonPython MODULE femriCommonPythonInit.cxx)
  TARGET_LINK_LIBRARIES(femriCommonPythonD femriCommon vtkCommonPythonD vtkCommon vtkImagingPythonD vtkImaging)
  TARGET_LINK_LIBRARIES (femriCommonPython femriCommonPythonD)
  IF(WIN32 AND NOT CYGWIN)
     SET_TARGET_PROPERTIES(femriCommonPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriCommonPythonD)
  INSTALL_TARGETS(${FEMRI_INSTALL_LIB_DIR} femriCommonPython)
ENDIF (FEMRI_WRAP_PYTHON)
