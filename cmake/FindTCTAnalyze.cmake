# - Find TCTAnalyze installation
# This module tries to find the TCTAnalyze source code on your system.
# Once done this will define
#
#  TCT_FOUND - system has TCTAnalyze installed

MESSAGE(STATUS "Looking for TCTAnalyze...")

FIND_PATH(TCT_INCLUDE
NAMES TCTScan.h MeasureWF.h
PATHS ${PROJECT_SOURCE_DIR}/external/TCTAnalyze.V1.0
)

FIND_PATH(TCT_LIBRARY
NAMES libTCTAnalyse.so
PATHS ${PROJECT_SOURCE_DIR}/external/TCTAnalyze.V1.0
)

IF (TCT_LIBRARY)
  IF(TCT_INCLUDE)
    set(TCT_FOUND TRUE)
    MESSAGE(STATUS "Found TCTAnalyze: ${TCT_INCLUDE}")
  ELSE(TCT_INCLUDE)
    set(TCT_FOUND FALSE)
    MESSAGE(FATAL_ERROR "TCTAnalyze headers NOT FOUND. Make sure to install the headers.")
  ENDIF(TCT_INCLUDE)
ELSE (TCT_LIBRARY)
    set(TCT_FOUND FALSE)
    MESSAGE(FATAL_ERROR "TCTAnalyze library NOT FOUND. Make sure to install the library.")
ENDIF (TCT_LIBRARY)

set(TCTAnalyze_INCLUDE
    ${TCT_INCLUDE}
)
set(TCTAnalyze_LIBRARIES
    "-L${TCT_LIBRARY} -lTCTAnalyse")

