# -*- coding: utf-8 -*-
import os.path


def get_filename(path: str):
    _, filename = os.path.split(path)
    return filename
