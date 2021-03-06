PROJECT(FEMRI_SCRIPTS)

SET(SCRIPTS_SRCS
  femri2dcylinderkspace.py
  femri2dboxkspace.py
  femri3dcylinderkspace.py
  femri3dboxkspace.py
  femri3dspherekspace.py
  femriemptykspace.py
  femrikspace.py
  femrikspacecenter.py
  femrikspaceerror.py
  femrikspacetoimage.py
  femrikspaceshift.py
  femrikspacezeropadding.py
  femriimagetokspace.py
  femrinumericmeshkspace.py
  femrinumericsurfacekspace.py
  femrisymbolicmeshkspace.py
  femrithicknessmeasurement.py
  )

SET(MODULE_SRCS
  femriscripts.py
  )

IF(NOT FEMRI_SCRIPTS_INSTALL_BIN_DIR)
  SET(FEMRI_SCRIPTS_INSTALL_BIN_DIR ${FEMRI_SCRIPTS_INSTALL_ROOT}/bin)
ENDIF(NOT FEMRI_SCRIPTS_INSTALL_BIN_DIR)

IF(NOT FEMRI_SCRIPTS_INSTALL_LIB_DIR)
  SET(FEMRI_SCRIPTS_INSTALL_LIB_DIR ${FEMRI_SCRIPTS_INSTALL_ROOT}/lib/femri)
ENDIF(NOT FEMRI_SCRIPTS_INSTALL_LIB_DIR)

FOREACH (SCRIPT_FILE  ${BASE_SRCS})
  CONFIGURE_FILE(${FEMRI_SCRIPTS_SOURCE_DIR}/${SCRIPT_FILE} ${FEMRI_SCRIPTS_BINARY_DIR}/${SCRIPT_FILE} COPYONLY)
ENDFOREACH (SCRIPT_FILE)

SET (STRIPPED_SCRIPTS_SRCS )
FOREACH (SCRIPT_FILE  ${SCRIPTS_SRCS})
  CONFIGURE_FILE(${FEMRI_SCRIPTS_SOURCE_DIR}/${SCRIPT_FILE} ${FEMRI_SCRIPTS_BINARY_DIR}/${SCRIPT_FILE} COPYONLY)
  STRING(REGEX REPLACE ".py" "" STRIPPED_SCRIPT_FILE ${SCRIPT_FILE})
  CONFIGURE_FILE(${FEMRI_SCRIPTS_SOURCE_DIR}/${SCRIPT_FILE} ${FEMRI_SCRIPTS_BINARY_DIR}/${STRIPPED_SCRIPT_FILE} COPYONLY)
  SET (STRIPPED_SCRIPTS_SRCS ${STRIPPED_SCRIPTS_SRCS} ${STRIPPED_SCRIPT_FILE})
ENDFOREACH (SCRIPT_FILE)

INSTALL_FILES(${FEMRI_SCRIPTS_INSTALL_LIB_DIR} .py ${SCRIPTS_SRCS} ${MODULE_SRCS})
INSTALL_PROGRAMS(${FEMRI_SCRIPTS_INSTALL_BIN_DIR} ${STRIPPED_SCRIPTS_SRCS})

