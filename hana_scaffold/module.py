# -*- coding: utf-8 -*-
import subprocess
from .enzyme import enzyme_validator
from .validator import file_exist_validator, file_list_exist_validator, integer_validator

_BIN_SEARCH_DIR = []


def __hana_exit(text):
    raise Exception(text)


def __hana_search_binary(bin_name: str):
    import platform
    import os

    def get_environ_paths():
        if platform.system() == 'Windows':
            return os.environ.get('PATH').split(';'), '.exe'
        return os.environ.get('PATH').split(':'), ''

    sys_path, bin_suffix = get_environ_paths()
    bin_filename = '{}{}'.format(bin_name, bin_suffix)
    search_dirs = [*_BIN_SEARCH_DIR, *sys_path]

    def is_executable(path: str):
        return os.path.isfile(path) and os.access(path, os.X_OK)

    for dir_path in search_dirs:
        expect_path = os.path.join(dir_path, bin_filename)
        if is_executable(expect_path):
            return expect_path
    return ''


def __hana_module(bin_name: str, params: dict, param_arg_map: dict):
    # Try to find the binary in environment.
    mod_bin_path = __hana_search_binary(bin_name)
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
    mod_proc = subprocess.Popen(mod_args)
    mod_proc.wait()


def add_search_path(path):
    _BIN_SEARCH_DIR.append(path)


def extract(fasta_path: str, mapping: list, output_prefix: str,
            enzyme: str, allele_table: str = None,
            enzyme_range: int = 1000, skip_range: bool = None,
            bam_mapq: int = 40, bam_skip_flag: bool = None,
            pairs_read_length: int = 170,
            search_buffer: int = 32, mapping_buffer: int = 512,
            threads: int = 1):
    __hana_module('hana_extract', locals(), {
        'fasta_path': ('-f', str, file_exist_validator),
        'mapping': ('-m', list, file_list_exist_validator),
        'allele_table': ('-a', str, file_exist_validator),
        'output_prefix': ('-o', str, None),
        'threads': ('-t', int, integer_validator),
        'enzyme': ('-e', str, enzyme_validator),
        'enzyme_range': ('-r', int, integer_validator),
        'bam_mapq': ('-q', int, integer_validator),
        'search_buffer': ('--search-buffer', int, integer_validator),
        'mapping_buffer': ('--mapping-buffer', int, integer_validator),
        'bam_skip_flag': ('--no-flag', bool, None),
        'skip_range': ('--no-range', bool, None)
    })
