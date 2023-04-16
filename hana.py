# -*- coding: utf-8 -*-
import os
import sys
import argparse
import automated.config as config


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('project_dir', metavar='PROJECT_DIR', help='Project directory')
    return parser.parse_args()


def main():
    args = parse_args()
    # Go to the directory and find the configure file, named in "hana_config.json".
    project_dir = os.path.abspath(args.project_dir)
    project_config = config.load_config(project_dir)
    # Check the environment.
    pass
    # Check resource of the config.
    project_files = config.check_config_input_files(project_config)

    return 0


if __name__ == '__main__':
    sys.exit(main())