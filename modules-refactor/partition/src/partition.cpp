#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "hmr_ui.hpp"

#include "partition.hpp"

bool vector_has_intersection(const HMR_CONTIG_ID_VEC& group_a, const HMR_CONTIG_ID_VEC& group_b)
{
    auto i = group_a.begin(), j = group_b.begin();
    while (i != group_a.end() && j != group_b.end())
    {
        if (*i < *j)
        {
            ++i;
        }
        else if (*i < *j)
        {
            ++j;
        }
        else
        {
            return true;
        }
    }
    return false;
}

bool partition_is_merge_valid(HMR_CONTIG_ID_VEC* group_a, HMR_CONTIG_ID_VEC* group_b, const HMR_ALLELE_TABLE& allele_table)
{
    HMR_CONTIG_ID_VEC* cluster_small, * cluster_large;
    if (group_a->size() < group_b->size())
    {
        cluster_small = group_a;
        cluster_large = group_b;
    }
    else
    {
        cluster_small = group_b;
        cluster_large = group_a;
    }
    //Go through the smaller cluster.
    for (int32_t contig_index : *cluster_small)
    {
        auto allele_record_iter = allele_table.find(contig_index);
        if (allele_record_iter == allele_table.end())
        {
            continue;
        }
        //Whether these two vectors contains the same items.
        if (vector_has_intersection(allele_record_iter->second, *cluster_large))
        {
            //When they has intersection, we cannot merge them
            return false;
        }
    }
    return true;
}

void partition_init_clusters(const HMR_NODES& nodes, const HMR_CONTIG_ID_VEC& invalid_nodes, CLUSTER_INFO& info)
{
    //Allocate the belongs array.
    info.belongs = static_cast<HMR_CONTIG_ID_VEC**>(malloc(sizeof(HMR_CONTIG_ID_VEC*) * nodes.size()));
    if (!info.belongs)
    {
        time_error(-1, "Failed to allocate memory for contig belong matrix.");
    }
    assert(info.belongs);
    size_t num_of_contigs = nodes.size(), invalid_pos = 0, num_of_invalids = invalid_nodes.size(), 
        num_of_clusters = num_of_contigs - num_of_invalids, cluster_pos = 0;
    info.clusters = static_cast<HMR_CONTIG_ID_VEC**>(malloc(sizeof(HMR_CONTIG_ID_VEC*) * num_of_clusters));
    if (!info.clusters)
    {
        time_error(-1, "Failed to allocate memory for clusters.");
    }
    assert(info.clusters);
    info.cluster_size = num_of_clusters;
    for (size_t i = 0; i < num_of_contigs; ++i)
    {
        //Check whether the node is marked as invalid.
        if (invalid_pos < num_of_invalids && i == invalid_nodes[invalid_pos])
        {
            info.belongs[i] = NULL;
            ++invalid_pos;
            continue;
        }
        //Create an cluster for the contig.
        HMR_CONTIG_ID_VEC* cluster = new HMR_CONTIG_ID_VEC();
        assert(cluster);
        cluster->push_back(static_cast<int32_t>(i));
        info.belongs[i] = cluster;
        info.clusters[cluster_pos] = cluster;
        ++cluster_pos;
    }
    assert(cluster_pos == num_of_clusters);
    //Prepare the link densities.
    info.link_densities.resize(num_of_contigs);
}

void partition_free_clusters(CLUSTER_INFO& info)
{
    if (info.merge)
    {
        delete[] info.merge;
    }
    for (int32_t i = 0; i < info.cluster_size; ++i)
    {
        delete info.clusters[i];
    }
    free(info.clusters);
    free(info.belongs);
}

void partition_edge_size_proc(uint64_t edge_size, void* user)
{
    CLUSTER_INFO* info = static_cast<CLUSTER_INFO*>(user);
    //Create the merge operation for each edge.
    info->merge = new CLUSTER_MERGE_OP[edge_size >> 1];
    if (!info->merge)
    {
        time_error(-1, "Failed to allocate density info memory.");
    }
    assert(info->merge);
    info->merge_size = 0;
}

void partition_edge_proc(HMR_EDGE_INFO* edges, int32_t edge_size, void* user)
{
    CLUSTER_INFO* info = static_cast<CLUSTER_INFO*>(user);
    //Loop through the partition, construct the edge and link densities.
    for (int32_t i = 0; i < edge_size; ++i)
    {
        const auto &edge = edges[i];
        //Build the edge densities.
        info->link_densities[edge.start].insert(std::make_pair(edge.end, edge.weights));
        //Check contig is invalid or not.
        if (info->belongs[edge.start] == NULL || info->belongs[edge.end] == NULL)
        {
            continue;
        }
        //Construct the merge operation.
        if (edge.start < edge.end)
        {
            CLUSTER_MERGE_OP& op = info->merge[info->merge_size];
            op.a = info->belongs[edge.start];
            op.b = info->belongs[edge.end];
            op.score = edge.weights;
            op.is_valid = true;
            ++info->merge_size;
        }
    }
}

inline double get_density(int32_t node_id, const CONTIG_LINK_DENSITY& node_density)
{
    const auto& iter = node_density.find(node_id);
    return iter == node_density.end() ? 0.0 : iter->second;
}

double get_total_linkage(HMR_CONTIG_ID_VEC* group_a, HMR_CONTIG_ID_VEC* group_b, const GRAPH_LINK_DENSITY& link_density)
{
    double total_linkage = 0.0;
    //Loop for all the nodes in the existed group.
    for (int32_t i : *group_a)
    {
        //Extract its node map.
        const auto& node_link_density = link_density[i];
        for (int32_t j : *group_b)
        {
            //Sum its weights.
            total_linkage += get_density(j, node_link_density);
        }
    }
    return total_linkage;
}

void partition_cluster(CLUSTER_INFO& cluster_info, int32_t num_of_groups)
{
    auto& merges = cluster_info.merge;
    auto& groups = cluster_info.clusters;
    //Loop until:
    //   - Nothing to merge
    //   - Cluster number reaches request.
    uint64_t op_counter = 0;
    size_t non_singleton_clusters = 0;
    const size_t non_skipped = (cluster_info.cluster_size >> 1);
    while (cluster_info.merge_size > 0 && cluster_info.cluster_size > static_cast<size_t>(num_of_groups))
    {
        //Find out the maximum weight in merge request.
        size_t max_weight_id = 0;
        {
            double max_score = merges[0].score;
            for (size_t i = 1; i < cluster_info.merge_size; ++i)
            {
                if (merges[i].score > max_score)
                {
                    max_weight_id = i;
                    max_score = merges[i].score;
                }
            }
        }
        //Get the top of the vector, which is the operation we are taking.
        CLUSTER_MERGE_OP& op = merges[max_weight_id];
        //We take this operation.
        op.is_valid = false;
        HMR_CONTIG_ID_VEC* group_b = op.b, * group_a = op.a;
        bool op_accept = true;
        //Check is this operation validate the allele table.
        if (cluster_info.allele_table && !partition_is_merge_valid(group_a, group_b, *cluster_info.allele_table))
        {
            //Mark the merge operation is not accepted.
            op_accept = false;
        }
        else
        {
            //What this magic?
            if (group_a->size() == 1)
            {
                ++non_singleton_clusters;
            }
            if (group_b->size() == 1)
            {
                ++non_singleton_clusters;
            }
            --non_singleton_clusters;
            //Change the belong pointers.
            for (const int32_t contig_id : *group_b)
            {
                cluster_info.belongs[contig_id] = group_a;
            }
            //Increase the group a, and sort it.
            group_a->insert(group_a->end(), group_b->begin(), group_b->end());
            std::sort(group_a->begin(), group_a->end());
            //Invalid all the merge operations with group a and b.
            for (size_t i = 0; i < cluster_info.merge_size; ++i)
            {
                auto* op_a = merges[i].a, * op_b = merges[i].b;
                //If the merge operation is related to group b, marked as invalid, will be removed.
                if (op_a == group_b || op_b == group_b || op_a == group_a || op_b == group_a)
                {
                    merges[i].is_valid = false;
                    continue;
                }
            }
            //Remove group b from clusters, and move group a to the end of the clusters.
            size_t group_offset = 0;
            for (size_t i = 0; i < cluster_info.cluster_size; ++i)
            {
                if (groups[i] == group_a || groups[i] == group_b)
                {
                    ++group_offset;
                    continue;
                }
                if (group_offset)
                {
                    groups[i - group_offset] = groups[i];
                }
            }
            //Decrease the cluster size.
            --cluster_info.cluster_size;
            //Recover group b.
            delete group_b;
            //Put group a at the end of the clusters.
            groups[cluster_info.cluster_size - 1] = group_a;
        }
        //Remove the merge operations which are not valid.
        size_t op_offset = 0;
        for (size_t i = 0; i < cluster_info.merge_size; ++i)
        {
            //When the merge operation is invalid, increase the offset, skip to next.
            if (!merges[i].is_valid)
            {
                ++op_offset;
                continue;
            }
            //Copy the current operation to several offset before.
            merges[i - op_offset] = merges[i];
        }
        //Update the merge operation size.
        cluster_info.merge_size -= op_offset;
        //Only generate new cluster merge request when the op is accepted.
        if (!op_accept)
        {
            continue;
        }
        //Create merge operations to new group a.
        //Calculate the cluster to new cluster offset.
        const double group_a_size = static_cast<double>(group_a->size());
        //printf("Added\n");
        for (size_t i = 0; i < cluster_info.cluster_size; ++i)
        {
            if (groups[i] == group_a)
            {
                continue;
            }
            //Calculate the weight of map.
            double average_linkage = get_total_linkage(groups[i], group_a, cluster_info.link_densities) / static_cast<double>(groups[i]->size()) / group_a_size;
            if (average_linkage <= 0.0)
            {
                continue;
            }
            //Save the current merge request.
            CLUSTER_MERGE_OP& op = merges[cluster_info.merge_size];
            op.a = groups[i];
            op.b = group_a;
            op.score = average_linkage;
            op.is_valid = true;
            ++cluster_info.merge_size;
        }
        //UI hints.
        ++op_counter;
        //Analyze the current clusters if enough merges occured.
        if (op_counter > non_skipped && non_singleton_clusters <= num_of_groups)
        {
            if (non_singleton_clusters == num_of_groups)
            {
                break;
            }
        }
        if (op_counter % 50 == 0)
        {
            time_print("%zu cluster(s) left.", cluster_info.cluster_size);
        }
    }
}

double partition_contig_cluster_linkage(int32_t contig_id, HMR_CONTIG_ID_VEC* cluster, const GRAPH_LINK_DENSITY& link_density, bool *has_linkage)
{
    //Get the contig id to cluster.
    double total_linkage = 0.0;
    //Extract its node map.
    const CONTIG_LINK_DENSITY& node_link_density = link_density[contig_id];
    for (int32_t group_id : *cluster)
    {
        auto group_id_iter = node_link_density.find(group_id);
        if (group_id_iter != node_link_density.end())
        {
            *has_linkage = true;
            total_linkage += group_id_iter->second;
        }
    }
    return total_linkage / static_cast<double>(cluster->size());
}

typedef struct RECOVER_LINKAGE
{
    HMR_CONTIG_ID_VEC* cluster;
    double average_linkage;
} RECOVER_LINKAGE;

bool partition_recover_linkage_greater(const RECOVER_LINKAGE& lhs, const RECOVER_LINKAGE& rhs)
{
    return lhs.average_linkage > rhs.average_linkage;
}

void partition_recover(const std::vector<HMR_CONTIG_ID_VEC*>& clusters, const HMR_CONTIG_ID_VEC& invalid_ids, const int32_t non_info_ratio, CLUSTER_INFO& info)
{
    RECOVER_LINKAGE *contig_linkages = new RECOVER_LINKAGE[clusters.size()];
    double non_info_ratio_f = static_cast<double>(non_info_ratio);
    assert(contig_linkages);
    //Loop and check all the cluster linkages.
    bool has_linkage;
    for (const int32_t contig_id : invalid_ids)
    {
        size_t linkage_size = 0;
        for (HMR_CONTIG_ID_VEC* cluster : clusters)
        {
            //Calculate the contig to cluster linkage.
            has_linkage = false;
            double linkage = partition_contig_cluster_linkage(contig_id, cluster, info.link_densities, &has_linkage);
            if (has_linkage)
            {
                contig_linkages[linkage_size] = RECOVER_LINKAGE{ cluster, linkage };
                ++linkage_size;
            }
        }
        if (linkage_size == 0)
        {
            continue;
        }
        //Calculate the pass ratio. (why it is a boolean equation?)
        std::sort(contig_linkages, contig_linkages + linkage_size, partition_recover_linkage_greater);
        const auto& best = contig_linkages[0];
        const auto& second_best = contig_linkages[1];
        if (best.average_linkage >= non_info_ratio_f &&
            ((linkage_size == 1) || (second_best.average_linkage == 0.0) || (best.average_linkage / second_best.average_linkage >= non_info_ratio_f)))
        {
            info.belongs[contig_id] = best.cluster;
        }
    }
    delete[] contig_linkages;
}
