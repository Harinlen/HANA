TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hana_draft

INCLUDEPATH = ../shared

SOURCES += \
    ../shared/hmr_algorithm.cpp \
    ../shared/hmr_args.cpp \
    ../shared/hmr_bin_file.cpp \
    ../shared/hmr_contig_graph.cpp \
    ../shared/hmr_path.cpp \
    ../shared/hmr_ui.cpp \
    src/args_draft.cpp \
    src/draft_mappings.cpp \
    src/main.cpp

HEADERS += \
    ../shared/hmr_global.hpp \
    ../shared/hmr_algorithm.hpp \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    ../shared/hmr_bin_file.hpp \
    ../shared/hmr_contig_graph.hpp \
    ../shared/hmr_contig_graph_type.hpp \
    ../shared/hmr_path.hpp \
    ../shared/hmr_ui.hpp \
    src/args_draft.hpp \
    src/draft_mappings.hpp \
    src/draft_mappings_type.hpp
