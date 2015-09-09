# - Find LCR installation
# This module tries to find the LCR source code on your system.
# Once done this will define
#
#  LCD_FOUND - system has LCR installed
#  LCR_INCLUDE_DIR - ~ the LCR include directory 
#  LCR_LIBRARY - libLCR

MESSAGE(STATUS "Looking for LCR...")

FIND_PATH(LCR_INCLUDE_DIR 
NAMES   LeCroy.h
PATHS   ${PROJECT_SOURCE_DIR}/external/LeCroyConverter/include
)

FIND_PATH(LCR_LIBRARY
NAMES libLeCroy.so
PATHS ${PROJECT_SOURCE_DIR}/external/LeCroyConverter/lib
)

IF (LCR_LIBRARY)
  IF(LCR_INCLUDE_DIR)
    set(LCR_FOUND TRUE)
    MESSAGE(STATUS "Found LCR: ${LCR_INCLUDE_DIR}, ${LCR_LIBRARY}")
  ELSE(LCR_INCLUDE_DIR)
    set(LCR_FOUND FALSE)
    MESSAGE(FATAL_ERROR "LCR headers NOT FOUND. Make sure to install the headers.")
  ENDIF(LCR_INCLUDE_DIR)
ELSE (LCR_LIBRARY)
    set(LCR_FOUND FALSE)
    MESSAGE(FATAL_ERROR "LCR library NOT FOUND. Make sure to install the library.")
ENDIF (LCR_LIBRARY)

set(LCR_INCLUDE_DIR
    ${LCR_INCLUDE_DIR}
)
set(LCR_LIBRARY
    "-L${LCR_LIBRARY} -lLeCroy")

