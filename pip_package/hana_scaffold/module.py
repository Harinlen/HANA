# -*- coding: utf-8 -*-
import os
from typing import List, Union
from .ext_binary import search_binary
from .breakpoint import run_command
from .enzyme import enzyme_validator
from .validator import file_exist_validator, file_list_exist_validator, integer_validator, float_validator


def __hana_module(bin_name: str, params: dict, param_arg_map: dict):
    # Try to find the binary in environment.
    mod_bin_path = search_binary(bin_name)
    if len(mod_bin_path) == 0:
        raise Exception('Failed to find HANA standard pipeline "{}", please check environment variable or reinstall HANA.'.format(bin_name))
    # Initial the mod arguments.
    mod_args = [mod_bin_path]
    # Loop in function parameters, construct the argument result.
    for param_name in params:
        if param_name in param_arg_map:
            arg, arg_type, validator = param_arg_map[param_name]
            # Check argument type.
            param_value = params[param_name]
            # Ignore the optional argument.
            if param_value is None:
                continue
            # Check argument type.
            if not isinstance(param_value, arg_type):
                raise Exception('Parameter "{}" is not type {}.'.format(param_name, str(arg_type)))
            # Validate the argument.
            if validator is not None:
                param_value = validator(param_value)
            # Set the argument.
            if isinstance(param_value, bool):
                if arg:
                    mod_args.append(arg)
            if isinstance(param_value, list):
                mod_args += [arg, *param_value]
            else:
                mod_args += [arg, str(param_value)]
    # Call the binary.
    run_command(mod_args)


def extract(contig_path: str, mapping: List[str], output_prefix: str,
            enzyme: Union[str, List[str]], allele_table: str = None,
            enzyme_range: int = 1000, skip_range: bool = None,
            bam_mapq: int = 40, bam_skip_flag: bool = None,
            pairs_read_length: int = 170,
            search_buffer: int = 32, mapping_buffer: int = 512,
            threads: int = 1):
    # Execute the binary.
    __hana_module('hana_extract', locals(), {
        'contig_path': ('-f', str, file_exist_validator),
        'mapping': ('-m', list, file_list_exist_validator),
        'allele_table': ('-a', str, file_exist_validator),
        'output_prefix': ('-o', str, None),
        'threads': ('-t', int, integer_validator),
        'enzyme': ('-e', Union[str, List[str]], enzyme_validator),
        'enzyme_range': ('-r', int, integer_validator),
        'bam_mapq': ('-q', int, integer_validator),
        'search_buffer': ('--search-buffer', int, integer_validator),
        'mapping_buffer': ('--mapping-buffer', int, integer_validator),
        'bam_skip_flag': ('--no-flag', bool, None),
        'skip_range': ('--no-range', bool, None)
    })
    # Check output file exists.
    nodes_path = '{}.hmr_nodes'.format(output_prefix)
    reads_path = '{}.hmr_reads'.format(output_prefix)
    file_exist_validator(nodes_path)
    file_exist_validator(reads_path)
    return nodes_path, reads_path


def draft(hana_nodes: str, hana_reads: str, output_prefix: str,
          allele_table: str = None, reads_buffer: int = 512,
          node_min_links: int = 3, node_min_re: int = 10,
          node_max_density: float = 2.0):
    __hana_module('hana_draft', locals(), {
        'hana_nodes': ('-n', str, file_exist_validator),
        'hana_reads': ('-r', str, file_exist_validator),
        'allele_table': ('-a', str, file_exist_validator),
        'output_prefix': ('-o', str, None),
        'reads_buffer': ('-b', int, integer_validator),
        'node_min_links': ('--min-links', int, integer_validator),
        'node_min_re': ('--min-re', int, integer_validator),
        'node_max_density': ('--max-link-density', float, float_validator),
    })
    # Check output file exists.
    edges_path = '{}.hmr_edges'.format(output_prefix)
    file_exist_validator(edges_path)
    return edges_path


def partition(hana_nodes: str, hana_edges: str, num_of_groups: int, output_prefix: str,
              hana_allele_table: str = None, edges_buffer: int = 512,
              node_non_informative_ratio: int = 3):
    __hana_module('hana_partition', locals(), {
        'hana_nodes': ('-n', str, file_exist_validator),
        'hana_edges': ('-e', str, file_exist_validator),
        'hana_allele_table': ('-a', str, file_exist_validator),
        'num_of_groups': ('-g', int, integer_validator),
        'output_prefix': ('-o', str, None),
        'edges_buffer': ('-b', int, integer_validator),
        'node_non_informative_ratio': ('--non-informative-ratio', int, integer_validator),
    })
    # This should generate a group of files.
    group_paths = ['{}_{}g{}.hmr_group'.format(output_prefix, ii+1, num_of_groups) for ii in range(num_of_groups)]
    file_list_exist_validator(group_paths)
    return group_paths


def ordering(hana_nodes: str, hana_edges: str, hana_group: str, output_path: str,
             threads: int = 1, edges_buffer: int = 512,
             ea_seed: int = 806, ea_mutation_rate: float = 0.2, ea_converge_gens: int = 5000,
             ea_max_gens: int = 1000000, ea_num_of_populations: int = 100):
    __hana_module('hana_ordering', locals(), {
        'hana_nodes': ('-n', str, file_exist_validator),
        'hana_edges': ('-e', str, file_exist_validator),
        'hana_group': ('-g', str, file_exist_validator),
        'output_path': ('-o', str, None),
        'threads': ('-t', int, integer_validator),
        'edges_buffer': ('-b', int, integer_validator),
        'ea_seed': ('-s', int, integer_validator),
        'ea_mutation_rate': ('--mutapb', float, float_validator),
        'ea_converge_gens': ('--ngen', int, integer_validator),
        'ea_max_gens': ('--max-gen', int, integer_validator),
        'ea_num_of_populations': ('--npop', int, integer_validator),
    })
    # Check whether output path is generated.
    file_exist_validator(output_path)
    return output_path


def orientation(hana_nodes: str, hana_reads: str, hana_seqs: List[str],
                reads_buffer: int = 512):
    __hana_module('hana_orientation', locals(), {
        'hana_nodes': ('-n', str, file_exist_validator),
        'hana_reads': ('-r', str, file_exist_validator),
        'hana_seqs': ('-s', list, file_list_exist_validator),
        'reads_buffer': ('-b', int, integer_validator),
    })
    # Check whether the orientation decided files exist.
    chromo_paths = []
    for seq_path in hana_seqs:
        chromo_path, _ = os.path.splitext(seq_path)
        chromo_path = '{}.hmr_chromo'.format(chromo_path)
        file_exist_validator(chromo_path)
        chromo_paths.append(chromo_path)
    return chromo_paths


def build(contig_path: str, hana_chromos: List[str], output_prefix: str):
    __hana_module('hana_build', locals(), {
        'contig_path': ('-f', str, file_exist_validator),
        'hana_chromos': ('-c', list, file_list_exist_validator),
        'output_prefix': ('-o', str, None),
    })
