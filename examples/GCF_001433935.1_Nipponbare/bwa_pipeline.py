# -*- coding: utf-8 -*-
import hana_scaffold
from hana_scaffold.pipeline import DiploidPipeline


hana_scaffold.add_search_path(r'/home/saki/Documents/hmr/HANA/bin/')
pipeline = DiploidPipeline(threads=16)
pipeline(contig_path='contigs_sim.fasta',
         lib_files=[('SRR6470741_1.fastq.gz', 'SRR6470741_2.fastq.gz')],
         enzyme='MboI', output_dir='/home/saki/Documents/hana-results/bare_auto/',
         num_of_groups=12)
