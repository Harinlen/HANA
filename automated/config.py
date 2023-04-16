# -*- coding: utf-8 -*-
import json
import os
import shutil
import platform

from .ui import time_exit, time_print

TEMPLATES = {
    'configure': 'hana_config.json',
    'status': 'hana_status.json'
}


def project_path(dir_path: str, template_id: str):
    return os.path.join(dir_path, TEMPLATES[template_id])


def load_json(dir_path: str, template_id: str):
    try:
        with open(project_path(dir_path, template_id), 'r', encoding='utf-8') as json_file:
            return json.load(json_file)
    except Exception:
        return {}


def save_json(content, dir_path: str, template_id: str):
    with open(project_path(dir_path, template_id), 'w', encoding='utf-8') as json_file:
        return json.dump(content, json_file)


def load_config(dir_path: str):
    return load_json(dir_path, 'configure')


def load_status(dir_path: str):
    return load_json(dir_path, 'status')


def dump_status(dir_path: str, status: dict):
    return save_json(status, dir_path, 'status')


def is_executable(path: str):
    return os.path.isfile(path) and os.access(path, os.X_OK)


def status_get_bwa(status: dict):
    if not is_executable(status['bwa_path']):
        return ''
    return status['bwa_path']


def status_get_samtools(status: dict):
    if not is_executable(status['samtools_path']):
        return ''
    return status['samtools_path']


def check_environment(config: dict, status: dict):
    if platform.system() == 'Windows':
        raise NotImplementedError
    else:
        search_dirs = os.environ.get('PATH').split(':')
        bin_suffix = ''
    # Load the path if config has paths settings.
    if 'paths' in config:
        # Load user settings to environment.
        path_config = config['paths']
        # Check whether it has directory settings.
        if 'dirs' in path_config:
            # Check directories as list of strings.
            env_paths = path_config['dirs']
            if not isinstance(env_paths, list):
                time_exit(1, '"dirs" should be a list of paths')
            for env_dir_path in env_paths:
                if not isinstance(env_dir_path, str):
                    time_exit(1, '"{}" is not a valid path.'.format(env_dir_path))
                if not os.path.isdir(env_dir_path):
                    time_print('Warning: "{}" is not a dir.'.format(env_dir_path))
            search_dirs = env_paths + search_dirs

    # Find binaries.
    def find_bin(bin_name: str):
        bin_name = '{}{}'.format(bin_name, bin_suffix)
        for dir_path in search_dirs:
            expect_path = os.path.join(dir_path, bin_name)
            if is_executable(expect_path):
                return expect_path
        return ''

    # Initial the path as empty.
    status['bwa_path'] = ''
    status['chromap_path'] = ''
    status['samtools_path'] = ''
    # Load the path if bwa or chromap path is manually set.
    if 'bwa' in config:
        status['bwa_path'] = config['bwa']
    if 'chromap' in config:
        status['chromap_path'] = config['chromap']
    if 'samtools' in config:
        status['samtools_path'] = config['samtools']
    # Find the bwa path and chromap path.
    if not is_executable(status['bwa_path']):
        status['bwa_path'] = find_bin('bwa')
    if not is_executable(status['chromap_path']):
        status['chromap_path'] = find_bin('chromap')
    if not is_executable(status['samtools_path']):
        status['samtools_path'] = find_bin('samtools')


# def check_config_input_files(config: dict):
#     reads_files = must_have(files, 'reads', 'Reads pair files not specify.')
#     if not isinstance(contig_file, str):
#         time_exit(1, 'Contig file path type does not match.')
#     if not os.path.isfile(contig_file):
#         time_exit(1, 'Contig file path is not a file.')
#     if not isinstance(reads_files, list):
#         time_exit(1, 'Reads pair list type does not match.')
#     for reads_file_path in reads_files:
#         if not os.path.isfile(reads_file_path):
#             time_exit(1, 'The following read file path is not a file:\n\t{}'.format(reads_file_path))
#     return files
#
