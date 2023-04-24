# -*- coding: utf-8 -*-
import sys
import subprocess
import importlib.util


def has_package(package_name: str):
    return importlib.util.find_spec(package_name) is not None


def pip_install(package_name: str):
    proc = subprocess.Popen([sys.executable, '-m', 'pip', 'install', package_name])
    proc.wait()


def ensure_package(package_name: str):
    if not has_package(package_name):
        pip_install(package_name)
