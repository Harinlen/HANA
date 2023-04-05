# HANA

High-performance Automated Next-generation AllHiC.

## Dependencies

- MSVC / GCC with C++ 11 support
- CMake
- Python 3.10

## Usage

### Arguments mode

```bash
hana -p PROJECT_DIR -f FASTA_PATH -1 [FORWARD_PAIRS] -2 [REVERSE_PAIRS] -g 12
```
Arguments explanation:
- `-p` or `--project-dir`: specify the project working directory.
- `-f` or `--fasta`: the config FASTA file. Its suffix should be either `.fasta` or `.fasta.gz`.
- `-1` or `--hic-forward`: the forward Hi-C FASTQ file. Its suffix should be `.fastq.gz`.
- `-2` or `--hic-reverse`: the reverse Hi-C FASTQ file. Its suffix should be `.fastq.gz`.
- `-g` or `--groups`: the number of chromosome groups.

Optional arguments:

- `-e` or `--enzyme`: the restriction sites used for the data. Can be a specific name or a small sequence. For more information, please check [HANA supported restriction sites](## Supported restriction sites). Default value is `HindIII`.
```bash
hana /path/to/directory
```

### Configuration file mode

## Supported restriction sites

