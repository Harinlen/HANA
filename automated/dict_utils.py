# -*- coding: utf-8 -*-
from .ui import time_exit


def must_have(source: dict, key: str, error_info: str):
    if key not in source:
        time_exit(error_info)
    return source[key]
