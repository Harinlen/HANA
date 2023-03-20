#ifndef DRAFT_MAP_TYPE
#define DRAFT_MAP_TYPE

#include <unordered_map>

#include "hmr_contig_graph_type.hpp"

#include "draft_fasta_type.hpp"

typedef std::unordered_map<std::string, int32_t> CONTIG_ID_MAP;

typedef struct MAPPING_COUNT
{
    int32_t pairs;              //Total pairs
    int32_t qualified_pairs;    //Valid pairs (>= quality)
} MAPPING_COUNT;
typedef std::unordered_map<uint64_t, MAPPING_COUNT> RAW_EDGE_MAP;

//Read position record.
typedef union
{
    struct {
        int32_t id;
        int32_t pos;
    } read;
    uint64_t data;
} READ_POS;
typedef std::unordered_set<uint64_t> READ_POS_SET;
typedef std::unordered_map<uint64_t, READ_POS_SET> READ_POS_MAP;

typedef struct
{
    const CONTIG_ID_MAP& contig_ids;
    const HMR_CONTIG_INVALID_SET& invalid_ids;
    const ENZYME_RANGES* contig_ranges;
    const uint8_t mapq;

    int32_t* contig_id_map; //Contig ID map: mapping the BAM contig to our own contig ID.
    int32_t contig_idx;     //Current processing Contig ID.
    FILE* reads_file;       //  BAM mappint file.
    char* output_buffer;    //  BAM mapping file buffer pointer.
    size_t output_offset,   //  BAM mapping file buffer offset.
        output_size;        //  BAM mapping file buffer total size.
    READ_POS_MAP records;   //Read records of the BAM files.
    RAW_EDGE_MAP edges;     //Edge raw records. (read pair counter)

} MAPPING_DRAFT_USER;

#endif // DRAFT_MAP_TYPE