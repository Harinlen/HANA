#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "partition.hpp"

CONTIG_ID_VECTORS partition_contigs(const HMR_CONTIGS &contigs, const HMR_EDGE_COUNTERS &edge_weights, int num_of_groups)
{
    //Prepare the contig sets.
    std::vector<CONTIG_ID_VECTOR *> contig_parts;
    contig_parts.reserve(contigs.size());
    std::unordered_set<int32_t> skipped_contigs;
    //Filter the count of RE.
    for(int32_t i=0, num_of_contigs=static_cast<int32_t>(contigs.size()); i<num_of_contigs; ++i)
    {
        if(contigs[i].enzyme_count < 25)
        {
            printf("%d\t%s\t%d\n", i, contigs[i].name, contigs[i].enzyme_count);
            skipped_contigs.insert(i);
        }
    }
    printf("%zu\n", skipped_contigs.size());
    //Assign contig set to NULLs for initialization.
    CONTIG_ID_VECTOR **contig_belongs = new CONTIG_ID_VECTOR *[contigs.size()];
    for(int32_t i=0, i_max = static_cast<int32_t>(contigs.size()); i<i_max; ++i)
    {
        //Create set for each contig.
        CONTIG_ID_VECTOR *contig_set = new CONTIG_ID_VECTOR();
        contig_set->push_back(i);
        //Insert the id_set to contig partitions.
        contig_parts.push_back(contig_set);
        contig_belongs[i] = contig_set;
    }
    //Loop until the contig parts reach the limitation.
    for(auto edge_info: edge_weights)
    {
        //Check the reference contigs has belong sets.
        int32_t contig_a = edge_info.edge.pos.start,
                contig_b = edge_info.edge.pos.end;
        if(skipped_contigs.find(contig_a) != skipped_contigs.end() ||
                skipped_contigs.find(contig_b) != skipped_contigs.end())
        {
            continue;
        }
        //Check the contig belongs is already assigned.
        CONTIG_ID_VECTOR *group_a = contig_belongs[contig_a],
                *group_b = contig_belongs[contig_b];
        if(group_a == group_b)
        {
            continue;
        }
        //Merge all the values in set b to set a.
        group_a->insert(group_a->end(), group_b->begin(), group_b->end());
        printf("Merge %d and %d -> %zu\n", contig_a, contig_b, group_a->size());
        for(int32_t id: *group_a)
        {
            printf("%d ", id);
        }
        printf("\n");
        for(int32_t contig_id: *group_b)
        {
            contig_belongs[contig_id] = group_a;
        }
        //Remove the group b from the set.
        contig_parts.erase(std::find(contig_parts.begin(), contig_parts.end(), group_b));
        delete group_b;
        if(static_cast<size_t>(num_of_groups) == contig_parts.size())
        {
            break;
        }
    }
    delete[] contig_belongs;
    //Create the contig result array and dump the data.
    CONTIG_ID_VECTORS result;
    result.reserve(contig_parts.size());
    for(CONTIG_ID_VECTOR *group_data: contig_parts)
    {
        //Dump the data to group vector.
        CONTIG_ID_VECTOR group;
        group.reserve(group_data->size());
        group.insert(group.end(), group_data->begin(), group_data->end());
        delete group_data;
        //Sort and save the group.
        std::sort(group.begin(), group.end());
        result.push_back(group);
    }
    return result;
}
