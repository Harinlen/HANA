TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hana_orientation

INCLUDEPATH = ../shared

SOURCES = \
    ../shared/hmr_args.cpp \
    ../shared/hmr_bin_file.cpp \
    ../shared/hmr_contig_graph.cpp \
    ../shared/hmr_path.cpp \
    ../shared/hmr_ui.cpp \
    src/args_orientation.cpp \
    src/main.cpp \
    src/orientation.cpp

HEADERS = \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    ../shared/hmr_bin_file.hpp \
    ../shared/hmr_contig_graph.hpp \
    ../shared/hmr_contig_graph_type.hpp \
    ../shared/hmr_path.hpp \
    ../shared/hmr_ui.hpp \
    src/args_orientation.hpp \
    src/orientation.hpp
