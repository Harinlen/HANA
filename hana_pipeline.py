# -*- coding: utf-8 -*-
import os
import sys
import argparse
import automated.config as config
import automated.ops as ops
from automated.ui import time_print


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('project_dir', metavar='PROJECT_DIR', help='Project directory')
    return parser.parse_args()


def main():
    args = parse_args()
    # Start up running the pipeline, print the version info.
    time_print("HANA Automated Pipeline v1.0")
    # Go to the directory and find the configure file, named in "hana_config.json".
    project_dir = os.path.abspath(args.project_dir)
    time_print("Running project at {}".format(project_dir))
    # Go to the project dir.
    os.chdir(project_dir)
    project_config = config.load_config(project_dir)
    # Check the environment.
    time_print("Loading environment...")
    project_status = config.load_status(project_dir)
    config.check_environment(project_config, project_status)
    # Run the operations.
    ops.run_pipeline(project_dir, project_config, project_status)
    return 0


if __name__ == '__main__':
    sys.exit(main())