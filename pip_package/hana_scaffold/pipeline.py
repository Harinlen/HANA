# -*- coding: utf-8 -*-
import os.path
from . import module as hana_module
from . import map as hana_map


class Pipeline:
    def __init__(self, **kwargs):
        self.settings = kwargs

    def __call__(self, *args, **kwargs):
        # Call the forward function.
        self.forward(*args, **kwargs)

    def forward(self, *args, **kwargs):
        pass


class DiploidPipeline(Pipeline):
    def __init__(self, **settings):
        super().__init__(**settings)

    def forward(self, contig_path: str, num_of_groups: int, output_dir: str, enzyme: str,
                lib_files: list = None, mapping_files: list = None, **kwargs):
        # Check library file paths or mapping files paths.
        if lib_files is None and mapping_files is None:
            raise Exception('Either Hi-C library files or mapped pairs files should be provided.')
        if lib_files is not None and mapping_files is not None:
            raise Exception('Only Hi-C library files or mapped pairs files need to be provided, not both.')
        # Based on the output directory, calculate the output prefix.
        dir_name = os.path.basename(output_dir)
        output_prefix = os.path.join(output_dir, dir_name)
        # Check to the target directory.
        os.chdir(output_dir)
        if lib_files is not None:
            # Map the Hi-C library, using BWA-MEM as default.
            mapper = hana_map.bwa_mem
            mapper_prefix = 'lib_{}.bwa_mem.bam'
            if 'mapper' in self.settings and self.settings['mapper'] == 'chromap':
                # Use chromap mapper.
                mapper = hana_map.chromap
                mapper_prefix = 'lib_{}.pairs'
            mapping_files = []
            for lib_id, (lib_forward_path, lib_reverse_path) in enumerate(lib_files):
                expect_lib_name = mapper_prefix.format(lib_id)
                if not os.path.isfile(expect_lib_name):
                    mapping_files.append(mapper(contig_path=contig_path,
                                                lib_forward_path=lib_forward_path,
                                                lib_reverse_path=lib_reverse_path,
                                                output_path=expect_lib_name,
                                                **self.settings))
                else:
                    mapping_files.append(expect_lib_name)
        # Call the pipelines.
        nodes_path, reads_path, _ = hana_module.extract(contig_path=contig_path,
                                                        mapping=mapping_files,
                                                        output_prefix=output_prefix, enzyme=enzyme,
                                                        **self.settings)
        edges_path = hana_module.draft(hana_nodes=nodes_path, hana_reads=reads_path, output_prefix=output_prefix,
                                       **self.settings)
        group_paths = hana_module.partition(hana_nodes=nodes_path, hana_edges=edges_path,
                                            output_prefix=output_prefix, **self.settings)
        seq_paths = []
        for group_path in group_paths:
            seq_path, _ = os.path.splitext(group_path)
            seq_path = '{}.hmr_seq'.format(seq_path)
            seq_paths.append(hana_module.ordering(hana_nodes=nodes_path, hana_edges=edges_path, hana_group=group_path,
                                                  output_path=seq_path, **self.settings))
        chromo_paths = hana_module.orientation(hana_nodes=nodes_path, hana_reads=reads_path, hana_seqs=seq_paths,
                                               **self.settings)
        hana_module.build(contig_path=contig_path, hana_chromos=chromo_paths, output_prefix=output_prefix)


class PolyploidPipeline(Pipeline):
    def __init__(self, **settings):
        super().__init__()
        self.pipeline = DiploidPipeline(**settings)

    def forward(self, contig_path: str, num_of_groups: int, allele_table: str, output_dir: str, enzyme: str,
                lib_files: list = None, mapping_files: list = None, **kwargs):
        self.pipeline.forward(contig_path, num_of_groups, output_dir, enzyme,
                              lib_files, mapping_files, allele_table=allele_table)

