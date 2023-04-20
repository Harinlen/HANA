# -*- coding: utf-8 -*-
from .ui import time_exit


def must_have(source: dict, key: str, error_info: str):
    if key not in source:
        time_exit(1, error_info)
    return source[key]


def get_var(source: dict, key: str, default: any = None):
    if key in source:
        return source[key]
    return default
