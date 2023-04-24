# -*- coding: utf-8 -*-
import json
import os
from automated.dict_utils import get_var, must_have
import platform

from .ui import time_exit, time_print

SCRIPT_DIR = os.path.dirname(__file__)
HANA_DIR = os.path.dirname(SCRIPT_DIR)
BIN_DIR = os.path.join(HANA_DIR, 'bin')


def get_environ_paths():
    if platform.system() == 'Windows':
        raise NotImplementedError
    return os.environ.get('PATH').split(':'), ''


def search_binary(tool_name: str, bin_suffix: str, search_dirs: list):
    bin_name = '{}{}'.format(tool_name, bin_suffix)
    for dir_path in search_dirs:
        expect_path = os.path.join(dir_path, bin_name)
        if is_executable(expect_path):
            return expect_path
    return ''


class Config:
    def __init__(self, dir_path: str):
        self.project_path = os.path.abspath(dir_path)
        _, self.project_dir = os.path.split(self.project_path)
        self.files = {}
        self.settings = {
            'threads': os.cpu_count(),
            'mapq': 40,
        }
        self.ops = {}
        self.__search_dirs, self.__bin_suffix = get_environ_paths()

    def get_hana_binary(self, name: str):
        bin_name = 'hana_{}{}'.format(name, self.__bin_suffix)
        bin_path = os.path.join(BIN_DIR, bin_name)
        if not is_executable(bin_path):
            time_exit(1, 'Failed to find HANA pipeline binary "{}"'.format(bin_name))
        return bin_path

    def get_files(self, key: str, must_exist: bool):
        if must_exist:
            return must_have(self.files, key, 'Failed to find "{}" files in config.'.format(key))
        return get_var(self.files, key, None)

    def get_global(self, key: str, must_exist: bool):
        if must_exist:
            return must_have(self.settings, key, 'Failed to find "{}" settings in global'.format(key))
        return get_var(self.settings, key, None)

    def parse_variable(self, var_key: str, variable: any, status: dict, must_exist: bool):
        if isinstance(variable, str):
            var_name = var_key if len(variable) == 1 else variable[1:]
            if variable.startswith('$'):
                # Load the variables from files.
                variable = self.get_files(var_name, must_exist)
            elif variable.startswith('!'):
                # Load the variables from global settings.
                variable = self.get_global(var_name, must_exist)
            elif variable.startswith('#'):
                # Get variable for project settings.
                var_name = var_name.lower()
                # Extract the variable from config.
                if not hasattr(self, var_name):
                    if must_exist:
                        time_exit(1, 'Failed to find "{}" variable in config'.format(var_name))
                    variable = None
                else:
                    variable = getattr(self, var_name)
            elif variable.startswith('-'):
                # Load the variables from status.
                if must_exist:
                    variable = must_have(status, var_name, 'Failed to find variable "{}" in status'.format(var_name))
                else:
                    variable = get_var(status, var_name)
            elif variable.startswith('[base]'):
                # Get the base path of the variable followed.
                variable = self.parse_variable('', variable[6:], status, True)
                variable, _ = os.path.splitext(variable)
            elif variable.startswith('[+]'):
                # We need to parse the following parts one by one as variable.
                variable = ''.join([self.parse_variable('', x.strip(), status, True) for x in variable[3:].split(',')])
        return variable

    def set_path(self, path_config: dict, status: dict):
        # Check binary paths.
        if not os.path.isdir(BIN_DIR):
            time_exit(1, 'Cannot find HANA binaries, please reinstall HANA.')

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
            self.__search_dirs = env_paths + self.__search_dirs

        def find_tool(tool_name: str):
            # Check whether we have a user preset path.
            if tool_name in path_config and is_executable(path_config[tool_name]):
                return path_config[tool_name]
            # Search the binary in system path.
            return search_binary(tool_name, self.__bin_suffix, self.__search_dirs)

        status['bwa_path'] = find_tool('bwa')
        status['chromap_path'] = find_tool('chromap')
        status['samtools_path'] = find_tool('samtools')


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
        return json.dump(content, json_file, indent='  ')


def load_config(dir_path: str, status: dict):
    config_dict = load_json(dir_path, 'configure')
    config = Config(dir_path)
    config.ops = get_var(config_dict, 'ops', {})
    config.files = get_var(config_dict, 'files', {})
    # Check files validation.
    if 'contigs' in config.files:
        contig_file = config.files['contigs']
        if not isinstance(contig_file, str):
            time_exit(1, 'Contig file path type does not match.')
        if not os.path.isfile(contig_file):
            time_exit(1, 'Contig file path is not a file.')
    if 'allele' in config.files:
        allele_file = config.files['allele']
        if not isinstance(allele_file, str):
            time_exit(1, 'Allele file path type does not match.')
        if not os.path.isfile(allele_file):
            time_exit(1, 'Allele file path is not a file.')

    def check_file_list(key: str, file_type: str):
        if key in config.files:
            key_files = config.files[key]
            if not isinstance(key_files, list):
                time_exit(1, '{} list type does not match.'.format(file_type))
            for key_file_path in key_files:
                if not isinstance(key_file_path, str):
                    time_exit(1, '"{}" is not a valid {} path.'.format(key_file_path, file_type))
                if not os.path.isfile(key_file_path):
                    time_exit(1, 'The following {} path is not a file:\n\t{}'.format(file_type, key_file_path))

    check_file_list('bams', 'BAM file')
    check_file_list('pairs', 'Pairs file')
    if 'mappings' in config.files:
        check_file_list('mappings', 'Mapping file')
    else:
        mappings_list = [*get_var(config.files, 'bams', []), *get_var(config.files, 'pairs', [])]
        config.files['mappings'] = mappings_list

    if 'reads' in config.files:
        reads_files = config.files['reads']
        if not isinstance(reads_files, list):
            time_exit(1, 'Reads pair list type does not match.')
        for reads_pair_paths in reads_files:
            if not isinstance(reads_pair_paths, list):
                time_exit(1, '"{}" is not a reads file pair list.'.format(reads_pair_paths))
            if len(reads_pair_paths) != 2:
                time_exit(1, '"{}" should only have 2 paths.'.format(reads_pair_paths))
            forward_path, reverse_path = reads_pair_paths
            if not os.path.isfile(forward_path):
                time_exit(1, 'The following read file path is not a file:\n\t{}'.format(forward_path))
            if not os.path.isfile(reverse_path):
                time_exit(1, 'The following read file path is not a file:\n\t{}'.format(reverse_path))

    # Update and check settings validation.
    config.settings.update(get_var(config_dict, 'global', {}))
    if not isinstance(config.settings['threads'], int):
        time_exit(1, '"globals" -> "threads" must be an integer')
    # Update and set paths for tools.
    config.set_path(get_var(config_dict, 'paths', {}), status)
    return config


def load_status(dir_path: str):
    return load_json(dir_path, 'status')


def dump_status(dir_path: str, status: dict):
    return save_json(status, dir_path, 'status')


def is_executable(path: str):
    return os.path.isfile(path) and os.access(path, os.X_OK)


def status_bin_path(status: dict, bin_name: str):
    path_key = '{}_path'.format(bin_name)
    if not is_executable(status[path_key]):
        return ''
    return status[path_key]
