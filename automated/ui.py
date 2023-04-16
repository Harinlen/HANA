# -*- coding: utf-8 -*-
import time


def time_print(*args, **kwargs):
    print("\033[32m{}\033[0m".format(time.strftime('[%H:%M:%S]', time.localtime(time.time()))), *args, **kwargs)


def time_exit(error_code: int, *args, **kwargs):
    time_print(*args, **kwargs)
    exit(error_code)
