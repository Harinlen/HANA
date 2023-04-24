# -*- coding: utf-8 -*-
import subprocess
import os
import sys
from automated.config import Config
from automated.ui import time_exit, time_print


class EnzymeDict(dict):
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            for key_item in key:
                super().__setitem__(key_item.upper(), value)
        else:
            super().__setitem__(key.upper(), value)


HANA_ENZYME_MAP = EnzymeDict()
HANA_ENZYME_MAP["TaqI", "Taq1"] = "TCGA"
HANA_ENZYME_MAP["HaeIII", "Hae3"] = "GGCC"
HANA_ENZYME_MAP["DpnI", "Dpn1", "DpnII", "Dpn2", "MboI", "Mbo1", "Sau3AI"] = "GATC"
HANA_ENZYME_MAP["AluI", "Alu1"] = "AGCT"
HANA_ENZYME_MAP["NlaIII", "Nla3", "FaeI", "Fae1", "FatI", "Fat1", "Hin1II", "Hsp92II"] = "CATG"
HANA_ENZYME_MAP["HpaII", "Hpa2"] = "CCGG"
HANA_ENZYME_MAP["FokI", "Fok1"] = "GGATG"
HANA_ENZYME_MAP["AaaI", "Aaa1"] = "CGGCG"
HANA_ENZYME_MAP["HgaI", "Hga1"] = "GACGC"
HANA_ENZYME_MAP["BglII", "Bgl2"] = "AGATCT"
HANA_ENZYME_MAP["EcoRV", "EcoR5"] = "GATATC"
HANA_ENZYME_MAP["EcoRI", "EcoR1"] = "GAATTC"
HANA_ENZYME_MAP["BamHI", "BamH1"] = "GGATCC"
HANA_ENZYME_MAP["HindIII", "Hind3"] = "AAGCTT"
HANA_ENZYME_MAP["KpnI", "Kpn1"] = "GGTACC"
HANA_ENZYME_MAP["XbaI", "Xba1"] = "TCTAGA"
HANA_ENZYME_MAP["XhoI", "Xho1"] = "CTCGAG"
HANA_ENZYME_MAP["SacI", "Sac1"] = "GAGCTC"
HANA_ENZYME_MAP["PstI", "Pst1"] = "CTGCAG"
HANA_ENZYME_MAP["SmaI", "Sma1"] = "CCCGGG"
HANA_ENZYME_MAP["PvuII", "Pvu2"] = "CAGCTG"
HANA_ENZYME_MAP["SalI", "Sal1"] = "GTCGAC"
HANA_ENZYME_MAP["ScaI", "Sca1"] = "AGTACT"
HANA_ENZYME_MAP["SpeI", "Spe1"] = "ACTAGT"
HANA_ENZYME_MAP["SphI", "Sph1"] = "GCATGC"
HANA_ENZYME_MAP["StuI", "Stu1"] = "AGGCCT"
HANA_ENZYME_MAP["NdeI", "Nde1"] = "CATATG"
HANA_ENZYME_MAP["NotI", "Not1"] = "GCGGCCGC"


def hana_enzyme_validator(enzyme: any):
    if not isinstance(enzyme, str):
        time_exit(1, '{} is not a valid restriction sites.'.format(enzyme))
    # Force enzyme to be upper.
    enzyme = enzyme.upper()
    # Check whether it is an alias name recorded in the map.
    if enzyme in HANA_ENZYME_MAP:
        enzyme = HANA_ENZYME_MAP[enzyme]
    # Check whether the enzyme is valid.
    for x in enzyme:
        if x not in ('A', 'C', 'T', 'G'):
            time_exit(1, 'Restriction site "{}" contains invalid base pair "{}"'.format(enzyme, x))
    return enzyme


def hana_standard_op(op_id: list, op: dict, config: Config, status: dict,
                     parameter_default: dict, parameter_map: dict, validator=None, post_result=None):
    # Find the HANA binaries.
    hana_command = op['command']
    time_print('Configuring HANA pipeline operation {}...'.format(hana_command))
    hana_bin_path = config.get_hana_binary(hana_command)
    # Loop and check operation config existed in parameter map.
    parameter = parameter_default.copy()
    parameter.update(op)
    del parameter['command']
    hana_args = {}
    for op_param in parameter:
        if op_param not in parameter_map:
            continue
        # Parse the parameter.
        hana_arg, hana_arg_type, hana_arg_validator, must_exist = parameter_map[op_param]
        op_var = config.parse_variable(op_param, parameter[op_param], status, must_exist)
        # Check variable and must have.
        if op_var is None and not must_exist:
            continue
        # Check whether the variable matches the requirements.
        if not isinstance(op_var, hana_arg_type):
            time_exit(1, '"{}" is not a {} type variable.'.format(op_var, hana_arg_type))
        # Validate the variable.
        if hana_arg_validator is not None:
            op_var = hana_arg_validator(op_var)
        # Check the variable type, construct the parameter.
        if isinstance(op_var, bool) or isinstance(op_var, int) or isinstance(op_var, float) or \
                isinstance(op_var, str) or isinstance(op_var, list):
            hana_args[hana_arg] = op_var
        else:
            hana_args[hana_arg] = str(op_var)
    if validator is not None:
        validator(hana_args)
    # Combine the hana arg dict into real arguments.
    cmd_args = [hana_bin_path]
    for arg_flag in hana_args:
        if isinstance(hana_args[arg_flag], bool):
            if hana_args[arg_flag]:
                cmd_args.append(arg_flag)
        elif isinstance(hana_args[arg_flag], list):
            cmd_args += [arg_flag, *hana_args[arg_flag]]
        elif isinstance(hana_args[arg_flag], str):
            cmd_args += [arg_flag, hana_args[arg_flag]]
        else:
            cmd_args += [arg_flag, str(hana_args[arg_flag])]
    # Call the HANA binaries.
    time_print("> {}".format(' '.join(cmd_args)))
    sys.stdout.flush()
    hana_proc = subprocess.Popen(cmd_args, stdout=sys.stdout, stderr=sys.stderr)
    hana_proc.wait()
    if hana_proc.returncode != 0:
        time_exit(1, 'Error happens when running HANA pipeline {}.'.format(hana_command))
    # When successful, check post result.
    if post_result is not None:
        post_result(hana_args)


def hana_op_extract(op_id: list, op: dict, config: Config, status: dict):
    def extract_validator(param: dict):
        if '-e' not in param:
            time_exit(1, 'Restriction site is not specify.')

    def extract_post_result(param: dict):
        # This would add multiple files to status.
        output_prefix = param['-o']

        def set_status_file(key: str, suffix: str):
            out_path = output_prefix + suffix
            if not os.path.isfile(out_path):
                time_exit(1, 'File {} not found.'.format(out_path))
            status[key] = out_path

        set_status_file('nodes', '.hmr_nodes')
        set_status_file('reads', '.hmr_reads')
        if '-a' in param:
            set_status_file('allele', '.hmr_allele_table')

    hana_standard_op(op_id, op, config, status, {
        "contigs": "$",
        "mappings": "$",
        "allele": "$",
        "output": "#PROJECT_DIR",
        "enzyme": "!",
        "threads": "!",
        "mapq": "!"
    }, {
        'contigs': ('-f', str, None, True),
        'mappings': ('-m', list, None, True),
        'allele': ('-a', str, None, False),
        'output': ('-o', str, None, True),
        'threads': ('-t', int, None, False),
        'enzyme': ('-e', str, hana_enzyme_validator, True),
        'enzyme_range': ('-r', int, None, False),
        'mapq': ('-q', int, None, False),
        'search-buffer': ('--search-buffer', int, None, False),
        'mapping-buffer': ('--mapping-buffer', int, None, False),
        'skip-check-flag': ('--no-flag', bool, None, False),
        'skip-check-range': ('--no-range', bool, None, False)
    }, extract_validator, extract_post_result)


def hana_op_draft(op_id: list, op: dict, config: Config, status: dict):
    def draft_post_result(param: dict):
        # This would add multiple files to status.
        output_prefix = param['-o']
        edge_path = output_prefix + '.hmr_edges'
        if not os.path.isfile(edge_path):
            time_exit(1, 'File {} not found.'.format(edge_path))
        status['edges'] = edge_path

    hana_standard_op(op_id, op, config, status, {
        "nodes": "-",
        "reads": "-",
        "output": "#PROJECT_DIR",
        "allele": "-",
    }, {
        'nodes': ('-n', str, None, True),
        'reads': ('-r', str, None, True),
        'allele': ('-a', str, None, False),
        'output': ('-o', str, None, True),
        'buffer-size': ('-b', int, None, False),
        'min-links': ('--min-links', int, None, False),
        'min-re': ('--min-re', int, None, False),
        'max-link-density': ('--max-link-density', float, None, False),
    }, post_result=draft_post_result)


def hana_op_partition(op_id: list, op: dict, config: Config, status: dict):
    def partition_post_result(param: dict):
        num_of_group = param['-g']
        output_prefix = param['-o']
        group_paths = []
        for ii in range(num_of_group):
            group_file_path = output_prefix + '_{}g{}.hmr_group'.format(ii + 1, num_of_group)
            if not os.path.isfile(group_file_path):
                time_exit(1, 'File {} not found.'.format(group_file_path))
            group_paths.append(group_file_path)
        status['groups'] = group_paths

    hana_standard_op(op_id, op, config, status, {
        "nodes": "-",
        "edges": "-",
        "output": "#PROJECT_DIR",
        "groups": "!",
    }, {
        'nodes': ('-n', str, None, True),
        'edges': ('-e', str, None, True),
        'allele': ('-a', str, None, False),
        'groups': ('-g', int, None, True),
        'output': ('-o', str, None, True),
        'buffer-size': ('-b', int, None, False),
        'non-informative-ratio': ('--non-informative-ratio', int, None, False),
    }, post_result=partition_post_result)


def hana_op_ordering(op_id: list, op: dict, config: Config, status: dict):
    def ordering_post_result(param: dict):
        output_path = param['-o']
        seq_paths = []
        if 'seq' in status:
            seq_paths = status['seq']
        if not isinstance(seq_paths, list):
            time_exit(1, "Status corrupted: 'seq' in status is not list.")
        if not os.path.isfile(output_path):
            time_exit(1, 'File {} not found.'.format(output_path))
        seq_paths.append(output_path)
        status['seq'] = seq_paths

    hana_standard_op(op_id, op, config, status, {
        "nodes": "-",
        "edges": "-",
        "group": "-",
        "output": "[+][base]-group,.hmr_seq",
        "threads": "!"
    }, {
        'nodes': ('-n', str, None, True),
        'edges': ('-e', str, None, True),
        'group': ('-g', str, None, True),
        'output': ('-o', str, None, True),
        'threads': ('-t', int, None, False),
        'buffer-size': ('-b', int, None, False),
        'seed': ('-s', int, None, False),
        'mutate-prob': ('--mutapb', float, None, False),
        'non-change-gen': ('--ngen', int, None, False),
        'max-gen': ('--max-gen', int, None, False),
        'num-of-pop': ('--npop', int, None, False),
    }, post_result=ordering_post_result)


def hana_op_orientation(op_id: list, op: dict, config: Config, status: dict):
    def orientation_post_result(param: dict):
        chromo_paths = []
        for seq_path in param['-s']:
            seq_base, _ = os.path.splitext(seq_path)
            chromo_path = seq_base+'.hmr_chromo'
            if not os.path.isfile(chromo_path):
                time_exit(1, 'File {} not found.'.format(chromo_path))
            chromo_paths.append(chromo_path)
        status['chromo'] = chromo_paths

    hana_standard_op(op_id, op, config, status, {
        "nodes": "-",
        "reads": "-",
        "seq": "-",
    }, {
        'nodes': ('-n', str, None, True),
        'reads': ('-r', str, None, True),
        'seq': ('-s', list, None, True),
        'buffer-size': ('-b', int, None, False),
    }, post_result=orientation_post_result)


def hana_op_build(op_id: list, op: dict, config: Config, status: dict):
    hana_standard_op(op_id, op, config, status, {
        "contigs": "$",
        "chromo": "-",
        "output": "#PROJECT_DIR",
    }, {
        'contigs': ('-f', str, None, True),
        'chromo': ('-c', list, None, True),
        'output': ('-o', str, None, True),
    })
