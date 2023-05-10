# -*- coding: utf-8 -*-
import os
import hana_scaffold
from hana_scaffold import module as scaffold
from hana_scaffold import breakpoint as bp
from hana_scaffold import map as mapping


hana_scaffold.add_search_path(r'/home/saki/Documents/hmr/HANA/bin/')
work_dir = '/home/saki/Documents/hana-results/bare_auto/'
contig_path = 'contigs_sim.fasta'
hic_lib_paths = [('SRR6470741_1.fastq.gz', 'SRR6470741_2.fastq.gz')]
restriction_enzyme = 'MboI'
work_prefix = 'bare_bam'
num_of_groups = 12
num_of_threads = 16

bp.init()
os.chdir(work_dir)
mapping_paths = []
for lib_id, (lib_forward_path, lib_reverse_path) in enumerate(hic_lib_paths):
    mapping_paths.append(mapping.chromap(contig_path=contig_path,
                                         lib_forward_path=lib_forward_path, lib_reverse_path=lib_reverse_path,
                                         pairs_path='lib_{}.pairs'.format(lib_id)))
nodes_path, reads_path, _ = scaffold.extract(contig_path=contig_path,
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
