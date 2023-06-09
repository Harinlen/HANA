# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import shutil
from pip_package.hana_scaffold.ext_binary import search_binary

HANA_DIR = os.path.abspath(os.path.dirname(__file__))
HANA_COMPONENTS = ['extract', 'draft', 'partition', 'ordering', 'orientation', 'build']
HANA_BIN_DIR = os.path.join(HANA_DIR, 'bin')
HANA_PY_DIR = os.path.join(HANA_DIR, 'automated')


def make_hana_unix_like(work_dir: str):
    work_dir = os.path.abspath(work_dir)
    print('Building HANA in {}'.format(work_dir))
    # Find the CMake directory.
    cmake_path = search_binary('cmake')
    if len(cmake_path) == 0:
        raise Exception('Failed to find "cmake", please make sure you can run "cmake" command.')
    make_path = search_binary('make')
    if len(make_path) == 0:
        raise Exception('Failed to find "make", please make sure you can run "make" command.')
    # Ensure the work dir exists.
    os.makedirs(work_dir, exist_ok=True)
    if not os.path.isdir(work_dir):
        raise Exception('Failed to create CMake work dir {}'.format(work_dir))
    # Call the cmake.
    os.chdir(work_dir)
    cmake_cmd = [cmake_path, os.path.join(HANA_DIR, 'modules')]
    print('> {}'.format(' '.join(cmake_cmd)))
    cmake_proc = subprocess.Popen(cmake_cmd)
    cmake_proc.wait()
    if cmake_proc.returncode != 0:
        print(1, 'Error happens during the cmake process.')
    # Do make inside the build directory.
    make_cmd = [make_path, '-j{}'.format(os.cpu_count())]
    print('> {}'.format(' '.join(make_cmd)))
    make_proc = subprocess.Popen(make_cmd)
    make_proc.wait()
    if make_proc.returncode != 0:
        raise Exception('Error happens during the make process.')
    # Make a directory named 'bin' under the source code directory.
    os.makedirs(HANA_BIN_DIR, exist_ok=True)
    if not os.path.isdir(HANA_BIN_DIR):
        raise Exception('Failed to create "bin" directory under HANA source code directory.')
    # Copy all the binaries to the bin dir.
    for hana_bin in HANA_COMPONENTS:
        shutil.copy2(os.path.join(work_dir, hana_bin, 'hana_{}'.format(hana_bin)),
                     os.path.join(HANA_BIN_DIR, 'hana_{}'.format(hana_bin)))


def main():
    # No matter how, make the HANA.
    make_hana_unix_like('cmake_build')
    # Check whether we have to install the binaries to a specific directory.
    if len(sys.argv) > 1:
        # Try to copy the binaries to the target directory.
        target_dir = sys.argv[1]
        if not os.path.isdir(target_dir):
            raise Exception('{} is not a directory'.format(target_dir))
        print('Installing hana...')
        for hana_bin in HANA_COMPONENTS:
            print('Copying hana_{}...'.format(hana_bin))
            shutil.copy2(os.path.join(HANA_BIN_DIR, 'hana_{}'.format(hana_bin)),
                         os.path.join(target_dir, 'hana_{}'.format(hana_bin)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
