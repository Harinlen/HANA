# -*- coding: utf-8 -*-
import os
import subprocess
import sys

from automated.ui import time_exit, time_print
from automated.config import Config, dump_status, status_bin_path
from automated.dict_utils import must_have, get_var
import automated.path_utils as path_uilts
import automated.hana_ops as hana_ops

KEY_LAST_OP = 'last_op'


def hana_op_map_mapper_args(op: dict, status: dict):
    # Check operation command type.
    mapper = 'chromap'
    if 'mapper' in op:
        if not isinstance(op['mapper'], str):
            time_exit(1, 'Invalid mapper {}'.format(op['mapper']))
        mapper = op['mapper']
    if mapper == 'chromap':
        raise NotImplementedError
    if mapper == 'bwa':
        bwa_path = status_bin_path(status, 'bwa')
        if len(bwa_path) == 0:
            time_exit(1, "No valid bwa executable found, please specify the path of bwa in config file.")
        samtools_path = status_bin_path(status, 'samtools')
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
        return mapper, (bwa_path, bwa_args, samtools_path)
    time_exit(1, 'Unknown mapper "{}" detected.'.format(mapper))


def hana_op_map_valid(op_id: list, op: dict, config: Config, status: dict):
    # Check whether the arguments are valid.
    hana_op_map_mapper_args(op, status)


def hana_op_map(op_id: list, op: dict, config: Config, status: dict):
    # Parse the contig file setting.
    contig_file = get_var(op, 'contigs', '')
    if not isinstance(contig_file, str):
        time_exit(1, 'Contig file path type does not match.')
    if len(contig_file) == 0:
        contig_file = config.get_files('contigs')
    if contig_file.startswith('$'):
        contig_file = config.get_files(contig_file[1:])
    if not os.path.isfile(contig_file):
        time_exit(1, 'Contig file path is not a file.')
    # Check whether the forward and reverse fastq files exist.
    reads_files = get_var(op, 'reads', [])
    if isinstance(reads_files, str) and reads_files.startswith('$'):
        reads_files = config.get_files(reads_files[1:])
    if len(reads_files) == 0:
        reads_files = config.get_files('reads')
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
    mapper, mapper_pack = hana_op_map_mapper_args(op, status)
    if mapper == 'chromap':
        raise NotImplementedError
    if mapper == 'bwa':
        bwa_path, bwa_args, samtools_path = mapper_pack

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
        output_prefix = path_uilts.get_filename(config.project_path)
        sam_view_args = ['view', '-hF', '256', '-']
        last_pair_id = get_var(status, 'last_pair_id', -1)
        for pair_id, (reads_forward, reads_reverse) in enumerate(reads_files):
            if pair_id < last_pair_id:
                continue
            bam_proc_args = [*bwa_args, contig_file, reads_forward, reads_reverse]
            output_bam_filename = '{}.map{}.pairs{}.bam'.format(output_prefix, op_id, pair_id)
            output_bam_path = os.path.join(config.project_path, output_bam_filename)
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
            dump_status(config.project_path, status)
    time_exit(1, 'Unknown mapper {}'.format(mapper))


def hana_op_custom(op_id: list, op: dict, config: Config, status: dict):
    pass


def hana_op_run_body(op_id_prefix: list, op_body: list, config: Config, status: dict):
    # Check whether the current operation is already complete.
    if len(op_id_prefix) < len(status[KEY_LAST_OP]):
        last_op_index = status[KEY_LAST_OP][len(op_id_prefix)]
    else:
        last_op_index = -1
    # Extend my prefix.
    op_id_prefix.append(-1)
    # Run all the operations.
    for op_id in range(last_op_index + 1, len(op_body)):
        op = op_body[op_id]
        op_id_prefix[-1] = op_id
        # Check current operation is valid or not.
        if not isinstance(op, dict):
            time_exit(1, 'Op {} detected in config file.'.format(op_id))
        # Check operation command.
        if 'command' not in op:
            time_exit(1, 'Op {} does not have a command.'.format(op_id))
        if 'enable' in op and op['enable'] is False:
            continue
        command = op['command']
        if command not in HANA_OPS:
            time_exit(1, 'Unknown op command {}'.format(command))
        # Call the operation.
        time_print("Running operation {}...".format(op_id))
        op_command, _ = HANA_OPS[command]
        op_command(op_id_prefix, op, config, status)
        # Update the current status.
        status[KEY_LAST_OP] = op_id_prefix
        dump_status(config.project_path, status)
    # After all the op is complete, remove the last value.
    op_id_prefix.pop()
    status[KEY_LAST_OP] = op_id_prefix
    dump_status(config.project_path, status)


def hana_op_loop(op_id: list, op: dict, config: Config, status: dict):
    # Extract the loop body from the op.
    loop_body = must_have(op, "body", 'Loop operation {} does not have body.'.format(op_id))
    # Extract the list to iterate.
    loop_list = must_have(op, "list", "Loop operation {} does not have iteration list.".format(op_id))
    loop_list = config.parse_variable('', loop_list, status, True)
    # Extract the list variable.
    loop_var = must_have(op, "var", "Loop operation {} does not have loop variable name.".format(op_id))
    time_print("Start loop in variable {}...".format(op['list']))
    # Start loop in the loop list.
    for ii, loop_var_item in enumerate(loop_list):
        # Set loop var in status.
        status[loop_var] = loop_var_item
        # Run loop body.
        hana_op_run_body(op_id, loop_body, config, status)


HANA_OPS = {
    'map': (hana_op_map, hana_op_map_valid),
    'extract': (hana_ops.hana_op_extract, None),
    'draft': (hana_ops.hana_op_draft, None),
    'partition': (hana_ops.hana_op_partition, None),
    'ordering': (hana_ops.hana_op_ordering, None),
    'orientation': (hana_ops.hana_op_orientation, None),
    'build': (hana_ops.hana_op_build, None),
    'loop': (hana_op_loop, None),
    'custom': (hana_op_custom, None),
}


def run_pipeline(config: Config, status: dict, project_status: dict):
    # Check last status or force restart.
    AWAYS_RESTART = True
    if KEY_LAST_OP not in status or AWAYS_RESTART:
        # status = project_status
        status[KEY_LAST_OP] = [-1]
    elif not isinstance(status[KEY_LAST_OP], list):
        time_print('Warning: Corrupted last operation breakpoint, restart the pipeline.')
        status = project_status
        status[KEY_LAST_OP] = [-1]
    # Check restart of pipeline.
    if len(status[KEY_LAST_OP]) == 0:
        time_print(end=" ")
        restart_result = input('Pipeline is already complete, restart? [y/n] ')
        if restart_result.upper() == 'Y':
            status = project_status
            status[KEY_LAST_OP] = [-1]
        else:
            time_exit(0, 'Pipeline already complete.')
    elif -1 < status[KEY_LAST_OP][0] < len(config.ops):
        time_print('Resume to last stopped step...')
    # Save the status right now.
    dump_status(config.project_path, status)
    # Run all the operations.
    hana_op_run_body([], config.ops, config, status)
    time_print('Pipeline execution completed.')
