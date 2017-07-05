#-------------------------------------------------
#
# Project created by QtCreator 2015-11-04T19:44:29
#
#-------------------------------------------------

QT       -= core gui

greaterThan(QT_MAJOR_VERSION, 4): QT -= widgets

TARGET = tbrowser
TEMPLATE = app

# For Windows
win32:LIBS += -LC:\root\root_v5.34.34\lib -llibCore -llibMathCore -llibTree -llibCint -llibRIO -llibNet -llibThread -llibGpad -llibGraf -llibRint -llibHist -llibMatrix
win32:INCLUDEPATH += C:\root\root_v5.34.34\include\

# For Linux
unix:LIBS += -L/home/maren/GEANT4/root_v5.34.34_install/lib/root -lCore -lMathCore -lTree -lCint -lRIO -lNet -lThread -lHist -lMatrix
unix:INCLUDEPATH += /home/maren/GEANT4/root_v5.34.34_install/include/root\


SOURCES += \
    src/tbrowser.cxx

