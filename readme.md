# HANA

High-performance Automated Next-generation AllHiC.

## Dependencies

- MSVC / GCC with C++ 11 support (MSVC >= 19.0, GCC >= 4.8.1)
- CMake (>=3.0)
- Python 3.10 or later

## Installation

1. Installed the compiler and CMake that meets the requirement of HANA.
2. Prepare the Python environment for HANA. Make sure when you run `python` it is calls a Python interpreter that meets the requirement.
3. `cd` to HANA source code directory and run `python hana_install.py`.

## Usage

HANA operates using a combination of project and configuration file mode. 

```bash
hana PROJECT_DIR
```

Optional arguments:

- `-e` or `--enzyme`: the restriction sites used for the data. Can be a specific name or a small sequence. For more information, please check [HANA supported restriction sites](## Supported restriction sites). Default value is `HindIII`.
```bash
hana /path/to/directory
```
