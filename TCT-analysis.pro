TARGET = tct-analysis
TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG += qt

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

DEFINES += USE_GUI=1

# For Windows
win32:LIBS += -LC:\root\root_v5.34.34\lib -llibCore -llibMathCore -llibTree -llibCint -llibRIO -llibNet -llibThread -llibGpad -llibGraf -llibRint -llibHist -llibMatrix
win32:INCLUDEPATH += C:\root\root_v5.34.34\include\

# For Linux
unix:QMAKE_CXXFLAGS += -std=c++11
unix:LIBS += -L/home/maren/GEANT4/root_v5.34.34_install/lib/root -lCore -lMathCore -lTree -lCint -lRIO -lNet -lThread -lHist -lMatrix
unix:INCLUDEPATH += /home/maren/GEANT4/root_v5.34.34_install/include/root\

# General
INCLUDEPATH += include

HEADERS += \
    include/modules/ModuleEdgeDepletion.h \
    include/modules/ModuleEdgeField.h \
    include/modules/ModuleEdgeFocus.h \
    include/modules/ModuleLaserAnalysis.h \
    include/modules/ModuleTopDepletion.h \
    include/modules/ModuleTopFocus.h \
    include/modules/ModuleTopMobility.h \
    include/acquisition.h \
    include/analysis.h \
    include/base.h \
    include/gui_consoleoutput.h \
    include/gui_folders.h \
    include/gui_sample.h \
    include/measurement.h \
    include/param.h \
    include/qdebugstream.h \
    include/sample.h \
    include/scanning.h \
    include/tct_config.h \
    include/TCTModule.h \
    include/TCTReader.h \
    include/util.h

SOURCES += \
    src/modules/ModuleEdgeDepletion.cc \
    src/modules/ModuleEdgeField.cc \
    src/modules/ModuleEdgeFocus.cc \
    src/modules/ModuleLaserAnalysis.cc \
    src/modules/ModuleTopDepletion.cc \
    src/modules/ModuleTopFocus.cc \
    src/modules/ModuleTopMobility.cc \
    src/acquisition.cc \
    src/analysis.cc \
    src/base.cc \
    src/gui_folders.cc \
    src/gui_sample.cc \
    src/main_gui.cxx \
    src/measurement.cc \
    src/sample.cc \
    src/scanning.cc \
    src/tct_config.cc \
    src/TCTModule.cc \
    src/TCTReader.cc \
    src/util.cc

FORMS += \
    forms/form_sample.ui \
    forms/form_parameters.ui \
    forms/form_folders.ui \
    forms/base.ui



