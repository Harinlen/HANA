#include <algorithm>
#include <cstdlib>
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
    auto a_iter = map.find(a);
    if(a_iter == map.end())
    {
        //Create a new map.
        CONTIG_LINK_DENSITY a_map;
        a_map.insert(std::make_pair(b, density));
        map.insert(std::make_pair(a, a_map));
    }
    else
    {
        a_iter->second.insert(std::make_pair(b, density));
    }
}

void partition_create_init_clusters(const HMR_CONTIGS &contigs, CLUSTER_INFO &cluster_info)
{
    //Find out the contigs
    cluster_info.cluster_size = contigs.size();
    cluster_info.clusters = new CONTIG_ID_VECTOR*[cluster_info.cluster_size];
    cluster_info.belongs = new CONTIG_ID_VECTOR*[cluster_info.cluster_size];
    for(int32_t i=0, i_max = static_cast<int32_t>(contigs.size()); i<i_max; ++i)
    {
        CONTIG_ID_VECTOR *cluster = new CONTIG_ID_VECTOR();
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
            time_print("Contig %s marked as Few REs (%d)", contigs[i].name, contigs[i].enzyme_count);
        }
    }
    return skipped_ids;
}

void partition_init_link_densities(const HMR_CONTIGS &contigs, const HMR_EDGE_MAP &edge_map,
                                   CLUSTER_INFO &cluster_info, const int32_t min_re, const int32_t max_link_density,
                                   const char *output_prefix)
{
    //Get the ignorance list of the few REs.
    time_print("Skipping contigs whose enzyme count is less than %d...", min_re);
    CONTIG_ID_SET skipped = partition_skip_few_res(contigs, min_re);
    //Prepare the total links.
    time_print("Skipping contigs likely from repetitive regions with multiplicity %d...", max_link_density);
    int32_t num_of_contigs = static_cast<int32_t>(contigs.size());
    //AllHiC uses weights as counts, which is not correct. It is actually weights.
    double total_links = 0.0;
    double *factors = new double[contigs.size()];
    for(int32_t i=0; i<num_of_contigs; ++i)
    {
        factors[i] = 0.0f;
    }
    for(const auto &edge_weight: edge_map)
    {
        HMR_EDGE edge;
        edge.data = edge_weight.first;
        //Calculate the sum of weights.
        double weight = edge_weight.second;
        total_links += weight;
        factors[edge.pos.start] += weight;
        factors[edge.pos.end] += weight;
    }
    //Calculate the average links, and calculate the contig factor.
    double avergage_links = total_links * 2.0 / static_cast<double>(num_of_contigs);
    for(int32_t i=0; i<num_of_contigs; ++i)
    {
        factors[i] /= avergage_links;
        if(factors[i] >= static_cast<double>(max_link_density))
        {
            skipped.insert(i);
            time_print("Contig %s marked as Repetitive (%.2lfx)", contigs[i].name, factors[i]);
        }
    }
    time_print("%zu contig(s) are skipped total.", skipped.size());
    //Update the link densities.
    time_print("Initiating merge requests and link densities...");
    cluster_info.merge = new MERGE_OP[edge_map.size()];
    cluster_info.merge_size = 0;
    cluster_info.link_densities.reserve(edge_map.size() << 1);
    for(const auto &link: edge_map)
    {
        //Update the link density.
        HMR_EDGE edge = direct_edge(link.first);
        if(hInSet(edge.pos.start, skipped) || hInSet(edge.pos.end, skipped))
        {
            continue;
        }
        //Calculate the density.
        double start_density = link.second / factors[edge.pos.start];
        insert_link_density(edge.pos.start, edge.pos.end, start_density, cluster_info.link_densities);
        insert_link_density(edge.pos.end, edge.pos.start, link.second / factors[edge.pos.end], cluster_info.link_densities);
        //Create the merge request.
        MERGE_OP &op = cluster_info.merge[cluster_info.merge_size];
        op.a = cluster_info.belongs[edge.pos.start];
        op.b = cluster_info.belongs[edge.pos.end];
        op.weight = start_density;
        op.is_valid = true;
        ++cluster_info.merge_size;
    }
    delete[] factors;
    time_print("%zu merge requests created.", cluster_info.merge_size);
    //Delete all the skipped sets.
    if(skipped.empty())
    {
        return;
    }
    char skipped_path[4096];
    sprintf(skipped_path, "%s_skipped.hmr_group", output_prefix);
    time_print("Saving %zu invalid contig(s) info to %s", skipped.size(), skipped_path);
    CONTIG_ID_VECTOR skipped_ids(skipped.begin(), skipped.end());
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
        if(i == skipped_ids[offset])
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
        const auto &node_link_density = link_density.find(i)->second;
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
        std::make_heap(merges, merges+cluster_info.merge_size, merge_op_comp);
        //Get the top of the vector, which is the operation we are taking.
        MERGE_OP &op = merges[0];
        //Merge all the contigs in op.b -> op.a
        op.is_valid = false;
        CONTIG_ID_VECTOR *group_b = op.b, *group_a = op.a;
        //Change the belongs of group b.
        for(const int32_t contig_id: *group_b)
        {
            cluster_info.belongs[contig_id] = group_a;
        }
        group_a->insert(group_a->end(), group_b->begin(), group_b->end());
        //We loop for all the rest of the merge operations.
        for(size_t i=1; i<cluster_info.merge_size; ++i)
        {
            CONTIG_ID_VECTOR *op_a = merges[i].a, *op_b = merges[i].b;
            //If the merge operation is related to group b, marked as invalid, will be removed.
            if(op_a == group_b || op_b == group_b || op_a == group_a || op_b == group_a)
            {
                merges[i].is_valid = false;
                continue;
            }
        }
        //Remove group b from clusters.
        remove_group_record(cluster_info, group_b);
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
