TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hana_extract

LIBS += -lz

INCLUDEPATH = ../shared

SOURCES += \
    ../shared/hmr_args.cpp \
    src/args_extract.cpp \
    src/main.cpp

HEADERS += \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    src/args_extract.hpp

