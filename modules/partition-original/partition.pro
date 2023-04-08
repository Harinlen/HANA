TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

TARGET = hana_partition

LIBS += -lz -lpthread

INCLUDEPATH = ../shared

HEADERS += \
    ../shared/hmr_allele.hpp \
    ../shared/hmr_allele_type.hpp \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    ../shared/hmr_bin_file.hpp \
    ../shared/hmr_contig_graph.hpp \
    ../shared/hmr_contig_graph_type.hpp \
    ../shared/hmr_global.hpp \
    ../shared/hmr_gz.hpp \
    ../shared/hmr_path.hpp \
    ../shared/hmr_text_file.hpp \
    ../shared/hmr_ui.hpp \
    src/args_partition.hpp \
    src/partition.hpp \
    src/partition_type.hpp

SOURCES += \
    ../shared/hmr_allele.cpp \
    ../shared/hmr_args.cpp \
    ../shared/hmr_bin_file.cpp \
    ../shared/hmr_contig_graph.cpp \
    ../shared/hmr_gz.cpp \
    ../shared/hmr_path.cpp \
    ../shared/hmr_text_file.cpp \
    ../shared/hmr_ui.cpp \
    src/args_partition.cpp \
    src/main.cpp \
    src/partition.cpp
