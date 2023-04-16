# -*- coding: utf-8 -*-
import json
import os.path
import shutil
from .ui import time_exit
from .dict_utils import must_have

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


def load_progress(dir_path: str):
    return load_json(dir_path, 'status')


def check_environment(config: dict):
    pass


def check_config_input_files(config: dict):
    # Check whether the config files exist.
    files = must_have(config, 'files', 'No input files specify.')
    contig_file = must_have(files, 'contigs', 'Contig file not specify.')
    reads_files = must_have(files, 'reads', 'Reads pair files not specify.')
    if not isinstance(contig_file, str):
        time_exit(1, 'Contig file path type does not match.')
    if not os.path.isfile(contig_file):
        time_exit(1, 'Contig file path is not a file.')
    if not isinstance(reads_files, list):
        time_exit(1, 'Reads pair list type does not match.')
    for reads_file_path in reads_files:
        if not os.path.isfile(reads_file_path):
            time_exit(1, 'The following read file path is not a file:\n\t{}'.format(reads_file_path))
    return files

