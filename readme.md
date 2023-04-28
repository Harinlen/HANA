# HANA

High-performance Automated Next-generation AllHiC.

## Dependencies

- MSVC / GCC with C++ 11 support (MSVC >= 19.0, GCC >= 4.8.1)
- CMake (>=3.0)
- Python 3.10 or later

## Installation

1. Installed the compiler and CMake that meets the requirement of HANA.
2. Prepare the Python environment for HANA. Make sure when you run `python` it is calls a Python interpreter that meets the requirement.
3. `cd` to HANA source code directory and run `python hana_install.py -i TARGET_DIR`, and make sure `TARGET_DIR` is in the environment path list.
4. Now `hana` command is ready for usage.

## Usage

HANA operates using a combination of project and configuration file mode. 

```bash
hana_scaffold PROJECT_DIR
```

In the `PROJECT_DIR` path, it should contain a file named `hana_config.json`. To compose a `hana_config.json`, please check the [HANA configuration file manual](config_manual.md).
