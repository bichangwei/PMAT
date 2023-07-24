#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Assembly of the error corrected data using Newbler.
"""

from log import Log
import shutil
import os
import subprocess
from check_file import remove_file, rename_file

log = Log()

def run_newbler(cpu, assembly_seq, output_path, mi=90, ml=40):

    # runAssembly_container = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../container/runAssembly.sif")

    if os.path.exists(os.path.abspath(os.path.dirname(__file__) + "/../container/runAssembly.sif")):
        runAssembly_container = os.path.join(os.path.abspath(os.path.dirname(__file__) + "/../container"), "runAssembly.sif")
    elif shutil.which('PMAT'):
        runAssembly_container = os.path.join(os.path.abspath(os.path.dirname(shutil.which('PMAT')) + "/../container"), "runAssembly.sif")
    else:
        log.Warning("runAssembly.sif installation error!")

    log.get_path(f'The path of Newbler : {runAssembly_container}')
    
    log.section_header("Reads assembly start...")

   # newbler_output = f'{output_path}/assembly_result'

    mount_output = os.path.join("/data", output_path.lstrip('/'))
    command = f'singularity exec -B {output_path}:{mount_output} -B {assembly_seq} {runAssembly_container} runAssembly -cpu {cpu} -het -sio -m -urt -large -s 100 -nobig -mi {mi} -ml {ml} -o {mount_output}/assembly_result {assembly_seq}'.split(' ')
    newbler_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Loop over stdout and stderr in real-time and print the output    
    # with open(f'{logfile}/newbler.log', 'w') as newblerlog:
    for line in iter(newbler_process.stdout.readline, b''):
        log.log(line.decode().strip().replace('\r', ''))
        # newblerlog.write(line.decode().strip().replace('\r', '')+'\n')

    for line in iter(newbler_process.stderr.readline, b''):
        log.log(line.decode().strip().replace('\r', ''))
        # newblerlog.write(line.decode().strip().replace('\r', '')+'\n')
    newbler_process.communicate()

    log.section_tail("Reads assembly end.")
    log.get_path(f'Assembly results path : {newbler_output}')

    rename_file(f'{newbler_output}/454AllContigs.fna', f'{newbler_output}/PMATAllContigs.fna')
    rename_file(f'{newbler_output}/454ContigGraph.txt', f'{newbler_output}/PMATContigGraph.txt')

    for rmopt in [os.path.join(newbler_output, opt) for opt in os.listdir(newbler_output) if opt.startswith('454')]:
        remove_file(rmopt)

    assembly_output = newbler_output

    return assembly_output

if __name__ == '__main__':
      run_newbler()

