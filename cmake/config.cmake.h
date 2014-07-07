#pragma once

// CMake uses config.cmake.h to generate config.h within the build folder.
#ifndef TCTANALYSIS_CONFIG_H
#define TCTANALYSIS_CONFIG_H

#define PACKAGE_NAME "@CMAKE_PROJECT_NAME@"
#define PACKAGE_VERSION "@TCTANALYSIS_LIB_VERSION@"
#define PACKAGE_STRING PACKAGE_NAME " " PACKAGE_VERSION

#define PACKAGE_SOURCE_DIR "@PROJECT_SOURCE_DIR@"

#endif
