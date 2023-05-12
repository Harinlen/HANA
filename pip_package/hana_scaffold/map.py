# -*- coding: utf-8 -*-
import os.path
from .ext_binary import search_binary
from .breakpoint import run_command, run_pipeline


def bwa_mem(contig_path: str, lib_forward_path: str, lib_reverse_path: str, output_path: str, threads: int, **kwargs):
    # Try to find the binary in environment.
    bwa_path = search_binary('bwa')
    if len(bwa_path) == 0:
        raise Exception('Failed to find BWA, please check environment variable or reinstall BWA.')
    samtools_path = search_binary('samtools')
    if len(samtools_path) == 0:
        raise Exception('Failed to find Samtools, please check environment variable or reinstall Samtools.')
    # Check whether the contig is already index.
    if not os.path.isfile(contig_path):
        raise FileNotFoundError('Contig file does not exist: {}'.format(contig_path))
    if not os.path.isfile(lib_forward_path):
        raise FileNotFoundError('Forward Hi-C library does not exist: {}'.format(lib_forward_path))
    if not os.path.isfile(lib_reverse_path):
        raise FileNotFoundError('Reverse Hi-C library does not exist: {}'.format(lib_reverse_path))

    # Check whether the contig FASTA file is already has index.
    def has_index():
        index_exts = ['amb', 'ann', 'bwt', 'pac', 'sa']
        for index_ext in index_exts:
            contig_index_path = '{}.{}'.format(contig_path, index_ext)
            if not os.path.isfile(contig_index_path):
                return False
        return True

    if not has_index():
        # Need to index the FASTA file first.
        run_command([bwa_path, 'index', contig_path])
    # Then execute the BWA MEM mapping.
    run_pipeline([[bwa_path, 'mem', '-SP5M', '-t', str(threads), contig_path, lib_forward_path, lib_reverse_path],
                  [samtools_path, 'view', '-hF', '256', '-'],
                  [samtools_path, 'sort', '-@', str(threads), '-o', output_path, '-T', 'tmp.ali']])
    if not os.path.isfile(output_path):
        raise FileNotFoundError('Mapping file not generated: {}'.format(output_path))
    return output_path


def chromap(contig_path: str, lib_forward_path: str, lib_reverse_path: str, output_path: str, **kwargs):
    # Try to find the binary in environment.
    chromap_path = search_binary('chromap')
    if len(chromap_path) == 0:
        raise Exception('Failed to find chromap, please check environment variable or reinstall chromap.')
    # Check whether the contig is already index.
    if not os.path.isfile(contig_path):
        raise FileNotFoundError('Contig file does not exist: {}'.format(contig_path))
    if not os.path.isfile(lib_forward_path):
        raise FileNotFoundError('Forward Hi-C library does not exist: {}'.format(lib_forward_path))
    if not os.path.isfile(lib_reverse_path):
        raise FileNotFoundError('Reverse Hi-C library does not exist: {}'.format(lib_reverse_path))

    # Check whether the contig FASTA file is already has index.
    fasta_index_path = '{}.chromap_index'.format(contig_path)
    if not os.path.isfile(fasta_index_path):
        run_command([chromap_path, '-i', '-r', contig_path, '-o', fasta_index_path])

    # Just run chromap.
    run_command([chromap_path, '--preset', 'hic', '-x', fasta_index_path, '-r', contig_path, '-1', lib_forward_path,
                 '-2', lib_reverse_path, '-o', output_path])
    if not os.path.isfile(output_path):
        raise FileNotFoundError('Mapping file not generated: {}'.format(output_path))
    return output_path
