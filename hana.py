# -*- coding: utf-8 -*-
import sys
import argparse
import automated.config as config


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--project-dir', help='Project directory')
    return parser.parse_args()


def main():
    args = parse_args()
    # Go to the directory and find the configure file, named in "hana_config.json".
    
    return 0


if __name__ == '__main__':
    sys.exit(main())