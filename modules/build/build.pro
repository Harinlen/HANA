TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hana_build

LIBS += -lz

INCLUDEPATH = ../shared

SOURCES += \
    ../shared/hmr_args.cpp \
    ../shared/hmr_bin_file.cpp \
    ../shared/hmr_bin_queue.cpp \
    ../shared/hmr_contig_graph.cpp \
    ../shared/hmr_fasta.cpp \
    ../shared/hmr_gz.cpp \
    ../shared/hmr_path.cpp \
    ../shared/hmr_seq.cpp \
    ../shared/hmr_text_file.cpp \
    ../shared/hmr_ui.cpp \
    src/args_build.cpp \
    src/main.cpp

HEADERS += \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    ../shared/hmr_bin_file.hpp \
    ../shared/hmr_bin_queue.hpp \
    ../shared/hmr_char.hpp \
    ../shared/hmr_contig_graph.hpp \
    ../shared/hmr_contig_graph_type.hpp \
    ../shared/hmr_fasta.hpp \
    ../shared/hmr_gz.hpp \
    ../shared/hmr_path.hpp \
    ../shared/hmr_seq.hpp \
    ../shared/hmr_text_file.hpp \
    ../shared/hmr_ui.hpp \
    src/args_build.hpp

