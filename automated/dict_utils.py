# -*- coding: utf-8 -*-
from .ui import time_exit


def must_have(source: dict, key: str, error_info: str):
    if key not in source:
        time_exit(error_info)
    return source[key]


def status_var(status: dict, key: str, default: any = None):
    if key in status:
        return status[key]
    return default
