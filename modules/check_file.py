#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of PMAT and is used to check file formats, 
file paths and to create, delete and modify file names.
"""

import os
from log import Log
import gzip

log = Log()

def check_file_format(filename, fmt=None):
    if os.path.splitext(filename)[1] == '.gz':
        raw_data =  gzip.open(filename, 'rt')
        firstline = raw_data.readline()
        raw_data.close()
        if firstline.startswith(">"):
            return "fasta"
        elif firstline.startswith("@"):
            return "fastq"
    else:
        raw_data = open(filename, 'r')
        firstline = raw_data.readline().strip()
        raw_data.close()
        if firstline.startswith(">"):
            return "fasta"
        elif firstline.startswith("@"):
            return "fastq"


def check_file_path(pathname, filename):
    if os.path.exists(os.path.join(pathname, filename)):
        log.Info(log.bold_dim(' >>> ') + log.purple(f'{filename} exists'))
    else:
        log.quit_with_error(f'{filename} not exists')

def mkdir_file_path(pathname):
    if os.path.exists(pathname):
        pass
    else:
        os.makedirs(pathname)
    
def remove_file(pathname):
    if os.path.exists(pathname):
        os.remove(pathname)
        
def rename_file(old_pathnm, new_pathnm=False):
    if os.path.exists(old_pathnm):
        os.rename(old_pathnm, new_pathnm)
