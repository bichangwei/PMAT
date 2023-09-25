#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Print log messages during the run.
"""

import sys
import datetime

class Log:
    def __init__(self):
        self.END = '\033[0m'
        self.BOLD = '\033[1m'
        self.RED = '\033[31m'
        self.GREEN = '\033[32m'
        self.PURPLE = '\033[35m'
        self.BLUE = '\033[34m'
        self.DIM = '\033[2m'

    def log(self, message=''):
        print(message, file=sys.stdout, flush=True, end='\n')

    def Info(self, message=''):
        time = self.get_timestamp()
        time_str = self.dim(time)
        print('[' + self.bold_blue(' INFO ') + time_str + ' ] ' + message, file=sys.stdout, flush=True, end='\n')

    def Error(self, message=''):
        time = self.get_timestamp()
        time_str = self.dim(time)
        print('[' + self.bold_purple(' ERROR ') + time_str + ' ] ' + message, file=sys.stderr, flush=True, end='\n')

    def Warning(self, message=''):
        time = self.get_timestamp()
        time_str = self.dim(time)
        print('[' + self.bold_purple(' WARNING ') + time_str + ' ] ' + message, file=sys.stderr, flush=True, end='\n')
        sys.exit()

    def section_header(self, text):
        time = self.get_timestamp()
        time_str = self.dim(time)
        header = self.bold_blue(text)
        print(header + '\n' + time_str, file=sys.stdout, flush=True)

    def section_tail(self, text):
        time = self.get_timestamp()
        time_str = self.dim(time)
        tail = self.bold_blue(text)
        print(tail + '\n' + time_str, file=sys.stdout, flush=True)


    def get_path(self, text):
        get_path = text
        print(get_path, file=sys.stdout, flush=True, end='\n')


    

    def bold(self, text):
        return self.BOLD + text + self.END

    def red(self, text):
        return self.RED + text + self.END

    def bold_red(self, text):
        return self.RED + self.BOLD + text + self.END

    def purple(self, text):
        return self.PURPLE + text + self.END

    def bold_purple(self, text):
        return self.BOLD + self.PURPLE + text + self.END

    def green(self, text):
        return self.GREEN + text + self.END

    def bold_green(self, text):
        return self.BOLD + self.GREEN + text + self.END

    def blue(self, text):
        return self.BLUE + text + self.END

    def bold_blue(self, text):
        return self.BLUE + self.BOLD + text + self.END

    def dim(self, text):
        return self.DIM + text + self.END

    def bold_dim(self, text):
        return self.BOLD + self.DIM + text + self.END

    def get_timestamp(self):
        return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())