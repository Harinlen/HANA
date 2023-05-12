# -*- coding: utf-8 -*-
import hana_scaffold
from hana_scaffold.pipeline import DiploidPipeline


hana_scaffold.add_search_path(r'/home/saki/Documents/hmr/HANA/bin/')
pipeline = DiploidPipeline(contig_path='/home/saki/Documents/dataset/bare/contigs_sim.fasta',
                           mapping_files=['/home/saki/Documents/bwa-results/bare/sample.bwa_mem.bam'],
                           enzyme='MboI', output_dir='/home/saki/Documents/hana-results/bare_bam_pip/',
                           groups=12, threads=16)
pipeline()
