TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lz -lpthread
INCLUDEPATH = ../shared

HEADERS += \
    ../shared/hmr_args.hpp \
    ../shared/hmr_args_types.hpp \
    ../shared/hmr_bam.hpp \
    ../shared/hmr_bgzf.hpp \
    ../shared/hmr_bin_file.hpp \
    ../shared/hmr_bin_queue.hpp \
    ../shared/hmr_char.hpp \
    ../shared/hmr_contig_graph.hpp \
    ../shared/hmr_contig_graph_type.hpp \
    ../shared/hmr_enzyme.hpp \
    ../shared/hmr_fasta.hpp \
    ../shared/hmr_global.hpp \
    ../shared/hmr_gz.hpp \
    ../shared/hmr_mapping.hpp \
    ../shared/hmr_mapping_type.hpp \
    ../shared/hmr_parallel.hpp \
    ../shared/hmr_path.hpp \
    ../shared/hmr_seq.hpp \
    ../shared/hmr_text_file.hpp \
    ../shared/hmr_thread_pool.hpp \
    ../shared/hmr_ui.hpp \
    src/args_draft.hpp \
    src/draft_fasta_index.hpp \
    src/draft_fasta_type.hpp \
    src/draft_map_type.hpp \
    src/draft_mapping.hpp

SOURCES += \
    ../shared/hmr_args.cpp \
    ../shared/hmr_bam.cpp \
    ../shared/hmr_bgzf.cpp \
    ../shared/hmr_bin_file.cpp \
    ../shared/hmr_bin_queue.cpp \
    ../shared/hmr_contig_graph.cpp \
    ../shared/hmr_enzyme.cpp \
    ../shared/hmr_fasta.cpp \
    ../shared/hmr_gz.cpp \
    ../shared/hmr_mapping.cpp \
    ../shared/hmr_path.cpp \
    ../shared/hmr_seq.cpp \
    ../shared/hmr_text_file.cpp \
    ../shared/hmr_ui.cpp \
    src/args_draft.cpp \
    src/draft_fasta_index.cpp \
    src/draft_mapping.cpp \
    src/main.cpp

