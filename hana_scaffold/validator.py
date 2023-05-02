# -*- coding: utf-8 -*-
import os
from typing import List


def file_exist_validator(path: str):
    if not isinstance(path, str):
        raise Exception('{} is not a valid path type'.format(type(path)))
    # File must exist.
    if not os.path.isfile(path):
        raise Exception('File does not exist: {}'.format(path))
    return path


def file_list_exist_validator(paths: List[str]):
    if not isinstance(paths, list):
        raise Exception('{} is not a path list.'.format(type(paths)))
    # Valid each item in paths.
    for filepath in paths:
        file_exist_validator(filepath)
    return paths


def integer_validator(value: int):
    if not isinstance(value, int):
        raise Exception('{} is not an integer'.format(value))
    return value
