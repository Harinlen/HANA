# -*- coding: utf-8 -*-
from . import _BIN_SEARCH_DIR


def search_binary(bin_name: str):
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
