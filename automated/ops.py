# -*- coding: utf-8 -*-
import os
import subprocess
import sys

from automated.ui import time_exit, time_print
from automated.config import dump_status, status_get_bwa, status_get_samtools
from automated.dict_utils import must_have, status_var
import automated.path_utils as path_uilts

KEY_LAST_OP = 'last_op'


def hana_op_map(op_id: int, op: dict, project_dir: str, config: dict, status: dict):
    # Check whether the config files exist.
    files = must_have(config, 'files', 'No input files specify.')
    # Extract contig files.
    contig_file = must_have(files, 'contigs', 'Contig file not specify.')
    if not isinstance(contig_file, str):
        time_exit(1, 'Contig file path type does not match.')
    if not os.path.isfile(contig_file):
        time_exit(1, 'Contig file path is not a file.')
    # Check whether the forward and reverse fastq files exist.
    reads_files = must_have(files, 'reads', 'Read file not specify.')
    if not isinstance(reads_files, list):
        time_exit(1, 'Reads pair list type does not match.')
    for reads_file_pair in reads_files:
        if not isinstance(reads_file_pair, list):
            time_exit(1, '{} should be a pair of file paths.'.format(reads_file_pair))
        if len(reads_file_pair) != 2:
            time_exit(1, '{} should only have two file paths.'.format(reads_file_pair))
        if not os.path.isfile(reads_file_pair[0]):
            time_exit(1, 'Reads file not found: {}'.format(reads_file_pair[0]))
        if not os.path.isfile(reads_file_pair[1]):
            time_exit(1, 'Reads file not found: {}'.format(reads_file_pair[1]))
    # Check operation command type.
    mapper = 'chromap'
    if 'mapper' in op:
        if not isinstance(op['mapper'], str):
            time_exit(1, 'Invalid mapper {}'.format(op['mapper']))
        mapper = op['mapper']
    if mapper == 'chromap':
        raise NotImplementedError
    if mapper == 'bwa':
        bwa_path = status_get_bwa(status)
        if len(bwa_path) == 0:
            time_exit(1, "No valid bwa executable found, please specify the path of bwa in config file.")
        samtools_path = status_get_samtools(status)
        if len(samtools_path) == 0:
            time_exit(1, "No valid samtools executable found, please specify the path of samtools in config file.")
        bwa_args = ['-SP5M', '-t', str(os.cpu_count())]
        if 'args' in op:
            user_args = op['args']
            if not isinstance(user_args, list):
                time_exit(1, 'BWA MEM arguments should be a list of string, but not {}'.format(user_args))
            if len(user_args) > 0:
                for x in user_args:
                    if not isinstance(x, str):
                        time_exit(1, '{} is not a valid parameter', x)
                bwa_args = user_args
        # Check contig file is indexed or not.
        def has_samtools_index():
            index_suffix = ['.amb', '.ann', '.pac', '.bwt', '.sa']
            for suffix in index_suffix:
                if not os.path.isfile(contig_file + suffix):
                    return False
            return True

        if not has_samtools_index():
            time_print('BWA indexing {}'.format(contig_file))
            sys.stdout.flush()
            bwa_index_proc = subprocess.Popen([bwa_path, 'index', contig_file], stderr=sys.stdout)
            bwa_index_proc.wait()
            if bwa_index_proc.returncode != 0:
                time_exit(1, 'Error happens when indexing {} with BWA.'.format(contig_file))
        bwa_args.insert(0, 'mem')
        output_prefix = path_uilts.get_filename(project_dir)
        sam_view_args = ['view', '-hF', '256', '-']
        last_pair_id = status_var(status, 'last_pair_id', -1)
        for pair_id, (reads_forward, reads_reverse) in enumerate(reads_files):
            if pair_id < last_pair_id:
                continue
            bam_proc_args = [*bwa_args, contig_file, reads_forward, reads_reverse]
            output_bam_filename = '{}.map{}.pairs{}.bam'.format(output_prefix, op_id, pair_id)
            output_bam_path = os.path.join(project_dir, output_bam_filename)
            sam_sort_args = ['sort', '-@', str(os.cpu_count()), '-o', output_bam_path, '-T', 'tmp.ali']
            # Construct the pipeline.
            time_print("Mapping {} and {} to {}".format(path_uilts.get_filename(reads_forward),
                                                        path_uilts.get_filename(reads_reverse),
                                                        path_uilts.get_filename(contig_file)))
            sys.stdout.flush()
            proc_bwa = subprocess.Popen([bwa_path, *bam_proc_args], stdout=subprocess.PIPE, stderr=sys.stdout)
            proc_sam_view = subprocess.Popen([samtools_path, *sam_view_args], stdin=proc_bwa.stdout, stdout=subprocess.PIPE, stderr=sys.stdout)
            proc_sam_sort = subprocess.Popen([samtools_path, *sam_sort_args], stdin=proc_sam_view.stdout, stdout=subprocess.PIPE, stderr=sys.stdout)
            proc_bwa.wait()
            proc_sam_sort.wait()
            sys.stdout.flush()
            if proc_bwa.returncode != 0:
                time_exit(1, 'Error happens when mapping reads file with BWA.')
            status['pair_id'] = pair_id
            dump_status(project_dir, status)
    time_exit(1, 'Unknown mapper {}'.format(mapper))


def hana_op_custom(op_id: int, op: dict, project_dir: str, config: dict, status: dict):
    pass


HANA_OP_MAP = {
    'map': hana_op_map,
    'custom': hana_op_custom,
}


def run_pipeline(dir_path: str, config: dict, status: dict):
    # Extract ops.
    if 'ops' not in config:
        time_exit(1, 'No operations found in config file.')
    ops = config['ops']
    # Check last status.
    if KEY_LAST_OP not in status:
        status[KEY_LAST_OP] = -1
    elif not isinstance(status[KEY_LAST_OP], int):
        time_print('Warning: Corrupted last operation breakpoint, restart the pipeline.')
        status[KEY_LAST_OP] = -1
    elif status[KEY_LAST_OP] >= len(ops):
        time_print(end=" ")
        result = input('Pipeline is already complete, restart? [y/n] ')
        if result.upper() == 'Y':
            status[KEY_LAST_OP] = -1
            dump_status(dir_path, status)
        else:
            time_exit(0, 'Pipeline already complete.')
    elif -1 < status[KEY_LAST_OP] < len(ops):
        time_print('Resume to last stopped step...')
    # Run all the operations.
    for op_id in range(status[KEY_LAST_OP] + 1, len(ops)):
        op = ops[op_id]
        # Skip the operation when we have last op.
        if op_id < status[KEY_LAST_OP]:
            continue
        # Check current operation is valid or not.
        if not isinstance(op, dict):
            time_exit(1, 'Op {} detected in config file.'.format(op_id))
        # Check operation command.
        if 'command' not in op:
            time_exit(1, 'Op {} does not have a command.'.format(op_id))
        command = op['command']
        if command not in HANA_OP_MAP:
            time_exit(1, 'Unknown op command {}'.format(command))
        # Call the operation.
        time_print("Running operation {}...".format(op_id))
        HANA_OP_MAP[command](op_id, op, dir_path, config, status)
        # Update the current status.
        status[KEY_LAST_OP] = op_id
        dump_status(dir_path, status)
    time_print('Pipeline execution completed.')
    status[KEY_LAST_OP] = len(ops)
    dump_status(dir_path, status)
