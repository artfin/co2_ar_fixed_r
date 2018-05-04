QT += core
QT -= gui

TARGET = co2_ar_fixed_r 
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -O2 -Wall -Wextra -g

INCLUDEPATH += -lgsl -lgscblas
INCLUDEPATH += "/usr/local/include/eigen3"

SOURCES += main.cpp \
    gear.cpp \
    vmblock.cpp \
    awp.cpp \
    basis.cpp \
    fgauss.cpp \
    trajectory.cpp \
    fgauss.cpp \
    matrix_euler.cpp \
    psp_pes.cpp \
    leg_arr.cpp \
    co2_ar_dipole.cpp

HEADERS += \
    gear.hpp \
    vmblock.hpp \
    basis.hpp \
    awp.hpp \
    constants.hpp \
    trajectory.hpp \
    tags.hpp \
    matrix_euler.hpp \
    psp_pes.hpp \
    leg_arr.hpp \
    co2_ar_dipole.hpp
