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
    def __init__(self, contig_path: str, groups: int, output_dir: str, enzyme: str,
                 lib_files: list = None, mapping_files: list = None, **settings):
        super().__init__(**settings)
        # Check library file paths or mapping files paths.
        if lib_files is None and mapping_files is None:
            raise Exception('Either Hi-C library files or mapped pairs files should be provided.')
        if lib_files is not None and mapping_files is not None:
            raise Exception('Only Hi-C library files or mapped pairs files need to be provided, not both.')
        # Save the parameters.
        self.contig_path = contig_path
        self.lib_files = lib_files
        self.mapping_files = mapping_files
        self.num_of_groups = groups
        self.enzyme = enzyme
        # Based on the output directory, calculate the output prefix.
        self.output_dir = output_dir
        dir_name = os.path.basename(self.output_dir)
        self.output_prefix = os.path.join(self.output_dir, dir_name)

    def forward(self):
        # Check to the target directory.
        os.chdir(self.output_dir)
        if self.lib_files is not None:
            # Map the Hi-C library, using BWA-MEM as default.
            mapper = hana_map.bwa_mem
            mapper_prefix = 'lib_{}.bwa_mem.bam'
            if 'mapper' in self.settings and self.settings['mapper'] == 'chromap':
                # Use chromap mapper.
                mapper = hana_map.chromap
                mapper_prefix = 'lib_{}.pairs'
            self.mapping_files = []
            for lib_id, (lib_forward_path, lib_reverse_path) in enumerate(self.lib_files):
                expect_lib_name = mapper_prefix.format(lib_id)
                if not os.path.isfile(expect_lib_name):
                    self.mapping_files.append(mapper(contig_path=self.contig_path,
                                                     lib_forward_path=lib_forward_path,
                                                     lib_reverse_path=lib_reverse_path,
                                                     output_path=expect_lib_name,
                                                     **self.settings))
                else:
                    self.mapping_files.append(expect_lib_name)
        # Call the pipelines.
        nodes_path, reads_path, _ = hana_module.extract(contig_path=self.contig_path,
                                                        mapping=self.mapping_files,
                                                        output_prefix=self.output_prefix, enzyme=self.enzyme,
                                                        **self.settings)
        edges_path = hana_module.draft(hana_nodes=nodes_path, hana_reads=reads_path, output_prefix=self.output_prefix,
                                       **self.settings)
        group_paths = hana_module.partition(hana_nodes=nodes_path, hana_edges=edges_path,
                                            output_prefix=self.output_prefix, **self.settings)
        seq_paths = []
        for group_path in group_paths:
            seq_path, _ = os.path.splitext(group_path)
            seq_path = '{}.hmr_seq'.format(seq_path)
            seq_paths.append(hana_module.ordering(hana_nodes=nodes_path, hana_edges=edges_path, hana_group=group_path,
                                                  output_path=seq_path, **self.settings))
        chromo_paths = hana_module.orientation(hana_nodes=nodes_path, hana_reads=reads_path, hana_seqs=seq_paths,
                                               **self.settings)
        hana_module.build(contig_path=self.contig_path, hana_chromos=chromo_paths, output_prefix=self.output_prefix)


class PolyploidPipeline(DiploidPipeline):
    def __init__(self, contig_path: str, groups: int, allele_table: str, output_dir: str, enzyme: str,
                 lib_files: list = None, mapping_files: list = None, **settings):
        super().__init__(contig_path, groups, output_dir, enzyme, lib_files, mapping_files, allele_table=allele_table,
                         **settings)
