#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cstdio>

#include "hmr_global.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_ui.hpp"

#include "partition.hpp"

inline HMR_EDGE direct_edge(uint64_t data)
{
    HMR_EDGE edge{};
    edge.data = data;
    return edge;
}

inline void insert_link_density(int32_t a, int32_t b, double density, MAP_LINK_DENSITY &map)
{
    map[a].insert(std::make_pair(b, density));
}

void partition_create_init_clusters(const HMR_CONTIGS &contigs, CLUSTER_INFO &cluster_info)
{
    //Find out the contigs
    cluster_info.cluster_size = contigs.size();
    cluster_info.clusters = new CONTIG_ID_VECTOR*[cluster_info.cluster_size];
    cluster_info.belongs = new CONTIG_ID_VECTOR*[cluster_info.cluster_size];
    if (!cluster_info.clusters || !cluster_info.belongs)
    {
        time_error(-1, "No enough memory to allocate cluster belonging memory.\n");
    }
    assert(cluster_info.clusters);
    assert(cluster_info.belongs);
    for(int32_t i=0, i_max = static_cast<int32_t>(cluster_info.cluster_size); i<i_max; ++i)
    {
        CONTIG_ID_VECTOR *cluster = new CONTIG_ID_VECTOR();
        assert(cluster);
        cluster->push_back(i);
        cluster_info.clusters[i] = cluster;
        cluster_info.belongs[i] = cluster;
    }
}

CONTIG_ID_SET partition_skip_few_res(const HMR_CONTIGS &contigs, const int32_t min_re)
{
    CONTIG_ID_SET skipped_ids;
    for(size_t i=0; i<contigs.size(); ++i)
    {
        if(contigs[i].enzyme_count < min_re)
        {
            skipped_ids.insert(static_cast<int32_t>(i));
        }
    }
    return skipped_ids;
}

void partition_skip_empty_links(const MAP_LINK_DENSITY& link_densities, CONTIG_ID_SET& skipped)
{
    for (size_t i = 0; i < link_densities.size(); ++i)
    {
        if (link_densities[i].empty())
        {
            skipped.insert(static_cast<int32_t>(i));
        }
    }
}

void partition_init_factors(int32_t num_of_contigs, LINK_FACTORS& factors_info)
{
    //Prepare the density info.
    factors_info.factors = new double[num_of_contigs];
    if (!factors_info.factors)
    {
        time_error(-1, "No enough memory for factor links");
    }
    //Reset the factors.
    for (int32_t i = 0; i < num_of_contigs; ++i)
    {
        factors_info.factors[i] = 0.0f;
    }
}

void partition_free_factors(LINK_FACTORS& factors_info)
{
    delete[] factors_info.factors;
}

void partition_factors_size_proc(uint64_t, void*)
{
}

void partition_factors_edge_proc(const HMR_EDGE_INFO* edge_info, void* user)
{
    if (edge_info->start == edge_info->end)
    {
        return;
    }
    LINK_FACTORS* factors_info = reinterpret_cast<LINK_FACTORS*>(user);
    //Check allele table (prune).
    if (factors_info->allele_mode)
    {
        //Check whether start is in the table.
        auto iter = factors_info->allele_table.find(edge_info->start);
        if (iter != factors_info->allele_table.end())
        {
            if (iter->second.find(edge_info->end) != iter->second.end())
            {
                return;
            }
        }
    }
    //Insert and build the factors.
    factors_info->avergage_links += edge_info->weights;
    factors_info->factors[edge_info->start] += edge_info->weights;
    factors_info->factors[edge_info->end] += edge_info->weights;
}

void partition_link_densities_size_proc(uint64_t edge_size, void* user)
{
    DENSITY_INFO* density_info = reinterpret_cast<DENSITY_INFO*>(user);
    //Prepare the merge request and link densities.
    density_info->cluster.merge = new MERGE_OP[edge_size];
    if (!density_info->cluster.merge)
    {
        time_error(-1, "Failed to allocate density info memory.");
    }
    assert(density_info->cluster.merge);
    density_info->cluster.merge_size = 0;
}

void partition_link_densities_edge_proc(const HMR_EDGE_INFO* edge_info, void* user)
{
    DENSITY_INFO* density_info = reinterpret_cast<DENSITY_INFO*>(user);
    //Check allele table (prune).
    if (density_info->link.allele_mode)
    {
        //Check whether start is in the table.
        auto iter = density_info->link.allele_table.find(edge_info->start);
        if (iter != density_info->link.allele_table.end())
        {
            if (iter->second.find(edge_info->end) != iter->second.end())
            {
                return;
            }
        }
    }
    CONTIG_ID_SET &skipped = density_info->link.skipped;
    double* factors = density_info->link.factors;
    //If any spot is skipped, ignore the weight.
    int32_t start_index = edge_info->start, end_index = edge_info->end;
    if (start_index == end_index || hInSet(start_index, skipped) || hInSet(end_index, skipped))
    {
        return;
    }
    //Calculate the density.
    double start_density = edge_info->weights / factors[start_index];
    //Build the contig link graph.
    insert_link_density(start_index, end_index, start_density, density_info->cluster.link_densities);
    insert_link_density(end_index, start_index, edge_info->weights / factors[end_index], density_info->cluster.link_densities);
    //Create the merge request.
    MERGE_OP& op = density_info->cluster.merge[density_info->cluster.merge_size];
    op.a = density_info->cluster.belongs[start_index];
    op.b = density_info->cluster.belongs[end_index];
    op.weight = start_density;
    op.is_valid = true;
    ++density_info->cluster.merge_size;
}

void partition_remove_skipped_contigs(LINK_FACTORS &contig_factors, CLUSTER_INFO & cluster_info, const char *output_prefix)
{
    char skipped_path[4097];
#ifdef _MSC_VER
    sprintf_s(skipped_path, 4096, "%s_skipped.hmr_group", output_prefix);
#else
    sprintf(skipped_path, "%s_skipped.hmr_group", output_prefix);
#endif
    time_print("Saving %zu invalid contig(s) info to %s", contig_factors.skipped.size(), skipped_path);
    CONTIG_ID_VECTOR skipped_ids(contig_factors.skipped.begin(), contig_factors.skipped.end());
    std::sort(skipped_ids.begin(), skipped_ids.end());
    hmr_graph_save_partition(skipped_path, skipped_ids);
    time_print("Removing skipped clusters...");
    //Reset the belongs pointer.
    for(int32_t i: skipped_ids)
    {
        delete cluster_info.belongs[i];
        cluster_info.belongs[i] = NULL;
    }
    //Delete the cluster reference from the vector.
    size_t offset = 0;
    for(int32_t i=0, i_max = static_cast<int32_t>(cluster_info.cluster_size); i<i_max; ++i)
    {
        //If the current cluster we have to skip, then skip it.
        if(offset < skipped_ids.size() && i == skipped_ids[offset])
        {
            ++offset;
            continue;
        }
        //Move the useful cluster.
        if(offset != 0)
        {
            cluster_info.clusters[i - offset] = cluster_info.clusters[i];
        }
    }
    cluster_info.cluster_size -= offset;
    time_print("%zu clusters have been removed.", offset);
}

bool merge_op_comp(const MERGE_OP &a, const MERGE_OP &b)
{
    return a.weight < b.weight;
}

inline double get_density(int32_t node_id, const CONTIG_LINK_DENSITY &node_density)
{
    const auto &iter = node_density.find(node_id);
    return iter == node_density.end() ? 0.0 : iter->second;
}

double get_total_linkage(CONTIG_ID_VECTOR *group_a, CONTIG_ID_VECTOR *group_b,
                       const MAP_LINK_DENSITY &link_density)
{
    double total_linkage = 0.0;
    //Loop for all the nodes in the existed group.
    for(int32_t i: *group_a)
    {
        //Extract its node map.
        const auto &node_link_density = link_density[i];
        for(int32_t j: *group_b)
        {
            //Sum its weights.
            total_linkage += get_density(j, node_link_density);
        }
    }
    return total_linkage;
}

void remove_group_record(CLUSTER_INFO &cluster_info, CONTIG_ID_VECTOR* cluster)
{
    //Remove the cluster from group records.
    auto &groups = cluster_info.clusters;
    for(size_t i=0; i<cluster_info.cluster_size; ++i)
    {
        if(groups[i] == cluster)
        {
            //Copy all the reset to the current position.
            for(size_t j=i+1; j<cluster_info.cluster_size; ++j)
            {
                groups[j-1] = groups[j];
            }
            --cluster_info.cluster_size;
            break;
        }
    }
    //Free the cluter.
    delete cluster;
}

bool is_operation_valid(CONTIG_ID_VECTOR* group_a, CONTIG_ID_VECTOR* group_b, const HMR_ALLELE_TABLE &allele_table)
{
    CONTIG_ID_VECTOR* less, * large;
    if (group_a->size() < group_b->size())
    {
        less = group_a;
        large = group_b;
    }
    else
    {
        less = group_b;
        large = group_a;
    }
    for (int32_t contig_index : *less)
    {
        auto contig_invalid_iter = allele_table.find(contig_index);
        if (contig_invalid_iter == allele_table.end())
        {
            continue;
        }
        //Loop for everything in large set.
        auto &contig_invalid_set = contig_invalid_iter->second;
        for (int32_t candidate_index : *large)
        {
            if (hInSet(candidate_index, contig_invalid_set))
            {
                return false;
            }
        }
    }
    //All test passed.
    return true;
}

void partition_cluster(CLUSTER_INFO &cluster_info, int32_t num_of_groups)
{
    auto &merges = cluster_info.merge;
    //Loop until:
    //   - Nothing to merge
    //   - Cluster number reaches request.
    int32_t op_counter = 0;
    while(cluster_info.merge_size > 0 &&
          cluster_info.cluster_size > static_cast<size_t>(num_of_groups))
    {
        //Find out the maximum weight in merge request.
        size_t max_weight_id = 0;
        {
            double max_weight = merges[0].weight;
            for (size_t i = 1; i < cluster_info.merge_size; ++i)
            {
                if (merges[i].weight > max_weight)
                {
                    max_weight_id = i;
                    max_weight = merges[i].weight;
                }
            }
        }
        //Get the top of the vector, which is the operation we are taking.
        MERGE_OP &op = merges[max_weight_id];
        //Merge all the contigs in op.b -> op.a
        op.is_valid = false;
        CONTIG_ID_VECTOR* group_b = op.b, * group_a = op.a;
        bool op_accept = true;
        //Check is this operation validate the allele table.
        if (cluster_info.allele_mode && !is_operation_valid(group_a, group_b, cluster_info.allele_table))
        {
            //Mark the op is not accept.
            op_accept = false;
        }
        else
        {
            //Change the belongs of group b.
            for (const int32_t contig_id : *group_b)
            {
                cluster_info.belongs[contig_id] = group_a;
            }
            group_a->insert(group_a->end(), group_b->begin(), group_b->end());
            //We loop for all the rest of the merge operations.
            for (size_t i = 0; i < cluster_info.merge_size; ++i)
            {
                CONTIG_ID_VECTOR* op_a = merges[i].a, * op_b = merges[i].b;
                //If the merge operation is related to group b, marked as invalid, will be removed.
                if (op_a == group_b || op_b == group_b || op_a == group_a || op_b == group_a)
                {
                    merges[i].is_valid = false;
                    continue;
                }
            }
            //Remove group b from clusters.
            remove_group_record(cluster_info, group_b);
        }
        //Remove the merge operations which are not valid.
        size_t op_offset = 0;
        for(size_t i=0; i<cluster_info.merge_size; ++i)
        {
            //When the merge operation is invalid, increase the offset, skip to next.
            if(!merges[i].is_valid)
            {
                ++op_offset;
                continue;
            }
            //Copy the current operation to several offset before.
            merges[i-op_offset] = merges[i];
        }
        //Update the merge operation size.
        cluster_info.merge_size -= op_offset;
        //Only generate new cluster merge request when the op is accepted.
        if (!op_accept)
        {
            continue;
        }
        //Calculate the cluster to new cluster offset.
        auto &groups = cluster_info.clusters;
        const double group_a_size = static_cast<double>(group_a->size());
        for(size_t i=0; i<cluster_info.cluster_size; ++i)
        {
            if(groups[i] == group_a)
            {
                continue;
            }
            //Calculate the weight of map.
            double average_linkage = get_total_linkage(groups[i], group_a, cluster_info.link_densities) /
                    static_cast<double>(groups[i]->size()) / group_a_size;
            if(average_linkage <= 0.0)
            {
                continue;
            }
            //Save the current merge request.
            MERGE_OP &op = merges[cluster_info.merge_size];
            op.a = groups[i];
            op.b = group_a;
            op.weight = average_linkage;
            op.is_valid = true;
            ++cluster_info.merge_size;
        }
        ++op_counter;
        if(op_counter == 50)
        {
            time_print("%zu cluster(s) left.", cluster_info.cluster_size);
            op_counter = 0;
        }
    }
}
