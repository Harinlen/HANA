# -*- coding: utf-8 -*-
import os
import hana_scaffold
from hana_scaffold import module as scaffold
from hana_scaffold import breakpoint as bp


hana_scaffold.add_search_path(r'/home/saki/Documents/hmr/HANA/bin/')
contig_path = r'/home/saki/Documents/dataset/bare/contigs_sim.fasta'
mapping_paths = [r'/home/saki/Documents/chromap-results/bare/h01.pairs']
restriction_enzyme = 'MboI'
work_prefix = r'/home/saki/Documents/hana-results/bare_pairs_pip/bare_pairs'
num_of_groups = 12
num_of_threads = 16

bp.init()
nodes_path, reads_path = scaffold.extract(contig_path=contig_path,
                                          mapping=mapping_paths,
                                          output_prefix=work_prefix, enzyme=restriction_enzyme, threads=num_of_threads)
edges_path = scaffold.draft(hana_nodes=nodes_path, hana_reads=reads_path, output_prefix=work_prefix)
group_paths = scaffold.partition(hana_nodes=nodes_path, hana_edges=edges_path, output_prefix=work_prefix,
                                 num_of_groups=num_of_groups)
seq_paths = []
for group_path in group_paths:
    seq_path, _ = os.path.splitext(group_path)
    seq_path = '{}.hmr_seq'.format(seq_path)
    seq_paths.append(scaffold.ordering(hana_nodes=nodes_path, hana_edges=edges_path, hana_group=group_path,
                                       output_path=seq_path, threads=num_of_threads))
chromo_paths = scaffold.orientation(hana_nodes=nodes_path, hana_reads=reads_path, hana_seqs=seq_paths)
scaffold.build(contig_path=contig_path, hana_chromos=chromo_paths, output_prefix='{}_build'.format(work_prefix))
bp.finish()
