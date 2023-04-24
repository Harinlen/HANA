# -*- coding: utf-8 -*-
import os
import stat
import subprocess
import sys
import shutil
from automated.ui import time_print, time_exit
from automated.package import ensure_package
from automated.config import get_environ_paths, search_binary

HANA_DIR = os.path.abspath(os.path.dirname(__file__))
HANA_COMPONENTS = ['extract', 'draft', 'partition', 'ordering', 'orientation', 'build']
HANA_BIN_DIR = os.path.join(HANA_DIR, 'bin')
HANA_PY_DIR = os.path.join(HANA_DIR, 'automated')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--work_dir', metavar='WORK_DIR', default='cmake_build',
                        help='CMake work directory')
    parser.add_argument('-i', '--install', metavar='INSTALL_DIR', default='', help='Install directory')
    return parser.parse_args()


def make_hana_bash(dir_path: str, hana_dir: str):
    hana_bash = os.path.join(dir_path, 'hana')
    try:
        with open(hana_bash, 'w') as hana_bash_file:
            hana_bash_file.write('\n'.join([
                '#!/bin/bash',
                '{} {} "$@"'.format(sys.executable, os.path.join(hana_dir, 'hana_pipeline.py'))
            ]))
    except IOError:
        time_exit(1, 'Failed to compose hana binary wrapper under {}'.format(dir_path))
    # Add executable.
    bash_stat = os.stat(hana_bash)
    os.chmod(hana_bash, bash_stat.st_mode | stat.S_IEXEC)


def make_hana_unix_like(work_dir: str):
    work_dir = os.path.abspath(work_dir)
    time_print('Building HANA in {}'.format(work_dir))
    # Find the CMake directory.
    paths, suffix = get_environ_paths()
    cmake_path = search_binary('cmake', suffix, paths)
    if len(cmake_path) == 0:
        time_exit(1, 'Failed to find "cmake", please make sure you can run "cmake" command.')
    make_path = search_binary('make', suffix, paths)
    if len(make_path) == 0:
        time_exit(1, 'Failed to find "make", please make sure you can run "make" command.')
    # Ensure the work dir exists.
    os.makedirs(work_dir, exist_ok=True)
    if not os.path.isdir(work_dir):
        time_exit(1, 'Failed to create CMake work dir {}'.format(work_dir))
    # Call the cmake.
    os.chdir(work_dir)
    cmake_cmd = [cmake_path, os.path.join(HANA_DIR, 'modules'), '.']
    time_print('> {}'.format(' '.join(cmake_cmd)))
    cmake_proc = subprocess.Popen(cmake_cmd)
    cmake_proc.wait()
    # Do make inside the build directory.
    make_cmd = [make_path, '-j{}'.format(os.cpu_count())]
    time_print('> {}'.format(' '.join(make_cmd)))
    make_proc = subprocess.Popen(make_cmd)
    make_proc.wait()
    # Make a directory named 'bin' under the source code directory.
    os.makedirs(HANA_BIN_DIR, exist_ok=True)
    if not os.path.isdir(HANA_BIN_DIR):
        time_exit(1, 'Failed to create "bin" directory under HANA source code directory.')
    # Copy all the binaries to the bin dir.
    for hana_bin in HANA_COMPONENTS:
        shutil.copy2(os.path.join(work_dir, hana_bin, 'hana_{}'.format(hana_bin)),
                     os.path.join(HANA_BIN_DIR, 'hana_{}'.format(hana_bin)))
    # Create the hana call bash.
    make_hana_bash(HANA_DIR, HANA_DIR)


def install_hana(target_dir: str):
    target_dir = os.path.abspath(target_dir)
    if not os.path.isdir(target_dir):
        time_exit(1, 'Install dir not exist: {}'.format(target_dir))
    time_print('Installing HANA to {}'.format(target_dir))
    hana_dir = os.path.join(target_dir, 'hana_pipeline')
    if os.path.isdir(hana_dir):
        time_print('Existed HANA installation detected, overwrite?', end='')
        response = input(' [y/n] ').lower()
        if response == 'n':
            return
        # Remove the directory.
        shutil.rmtree(hana_dir)
    # Create HANA directory under target directory.
    os.makedirs(hana_dir, exist_ok=True)
    if not os.path.isdir(hana_dir):
        time_exit(1, 'Failed to create HANA directory under {}'.format(target_dir))
    # Copy the automated directory and hana_pipeline.py to the target dir.
    shutil.copytree(HANA_PY_DIR, os.path.join(hana_dir, 'automated'))
    shutil.copy2(os.path.join(HANA_DIR, 'hana_pipeline.py'), os.path.join(hana_dir, 'hana_pipeline.py'))
    # Copy the binary directory to the target.
    shutil.copytree(HANA_BIN_DIR, os.path.join(hana_dir, 'bin'))
    # Construct a bash script to target dir.
    make_hana_bash(target_dir, hana_dir)


def main():
    # Check the dependencies.
    ensure_package('argparse')
    # Find the CMake binaries.
    args = parse_args()
    # No matter how, make the HANA.
    make_hana_unix_like(args.work_dir)
    # Ensure the working directory.
    if len(args.install) > 0:
        # Check whether the build is complete.
        install_hana(args.install)
    return 0


if __name__ == '__main__':
    sys.exit(main())
