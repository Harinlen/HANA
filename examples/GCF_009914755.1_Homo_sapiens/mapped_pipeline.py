# -*- coding: utf-8 -*-
import hana_scaffold
from hana_scaffold.pipeline import DiploidPipeline


hana_scaffold.add_search_path(r'/home/saki/Documents/hmr/HANA/bin/')
pipeline = DiploidPipeline(threads=16)
pipeline(contig_path='/home/saki/Documents/bwa-results/human/contigs_sim.fasta',
         mapping_files=[r'/home/saki/Documents/bwa-results/human/SRR13849427.bwa_mem.bam',
                        r'/home/saki/Documents/bwa-results/human/SRR13849428.bwa_mem.bam',
                        r'/home/saki/Documents/bwa-results/human/SRR13849429.bwa_mem.bam',
                        r'/home/saki/Documents/bwa-results/human/SRR13849430.bwa_mem.bam'],
         enzyme=['GATC', 'GAAT', 'GACT', 'GAGT', 'GATT'],
         weight_enzyme=['GATCGATC', 'GANTGATC', 'GANTANTC', 'GATCANTC'],
         output_dir='/home/saki/Documents/hana-results/human_pip/',
         num_of_groups=24)
