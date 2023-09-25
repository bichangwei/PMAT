#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Assembly of the error corrected data using runAssembly.
"""

from log import Log
import shutil
import os
import subprocess
from check_file import remove_file, rename_file, remove_dir

log = Log()

def run_Assembly(cpu, assembly_seq, output_path, mi=90, ml=40, mem=None):

    # runAssembly_container = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../container/runAssembly.sif")

    if os.path.exists(os.path.abspath(os.path.dirname(__file__) + "/../container/runAssembly.sif")):
        runAssembly_container = os.path.join(os.path.abspath(os.path.dirname(__file__) + "/../container"), "runAssembly.sif")
    elif shutil.which('PMAT'):
        runAssembly_container = os.path.join(os.path.abspath(os.path.dirname(shutil.which('PMAT')) + "/../container"), "runAssembly.sif")
    else:
        log.Warning("runAssembly.sif installation error!")

    log.get_path(f'The path of runAssembly : {runAssembly_container}')
    
    log.section_header("Reads assembly start...")


    mount_output = os.path.join("/data", output_path.lstrip('/'))
    runAssembly_output = f'{output_path}/assembly_result'
    if mem:
        command = f'singularity exec -B {output_path}:{mount_output} -B {assembly_seq} {runAssembly_container} runAssembly -cpu {cpu} -het -sio -m -urt -large -s 100 -nobig -mi {mi} -ml {ml} -o {mount_output}/assembly_result {assembly_seq}'.split(' ')
    else:
        command = f'singularity exec -B {output_path}:{mount_output} -B {assembly_seq} {runAssembly_container} runAssembly -cpu {cpu} -het -sio -urt -large -s 100 -nobig -mi {mi} -ml {ml} -o {mount_output}/assembly_result {assembly_seq}'.split(' ')
    runAssembly_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Loop over stdout and stderr in real-time and print the output    
    # with open(f'{logfile}/runAssembly.log', 'w') as runAssemblylog:
    for line in iter(runAssembly_process.stdout.readline, b''):
        log.log(line.decode().strip().replace('\r', ''))
        # runAssemblylog.write(line.decode().strip().replace('\r', '')+'\n')

    for line in iter(runAssembly_process.stderr.readline, b''):
        log.log(line.decode().strip().replace('\r', ''))
        # runAssemblylog.write(line.decode().strip().replace('\r', '')+'\n')
    runAssembly_process.communicate()

    log.section_tail("Reads assembly end.")
    log.get_path(f'Assembly results path : {runAssembly_output}')

    if os.path.exists(f'{runAssembly_output}/454AllContigs.fna') and os.path.exists(f'{runAssembly_output}/454ContigGraph.txt'):
        rename_file(f'{runAssembly_output}/454AllContigs.fna', f'{runAssembly_output}/PMATAllContigs.fna')
        rename_file(f'{runAssembly_output}/454ContigGraph.txt', f'{runAssembly_output}/PMATContigGraph.txt')

        for rmopt in [os.path.join(runAssembly_output, opt) for opt in os.listdir(runAssembly_output) if opt.startswith('454')]:
            remove_file(rmopt)
        remove_dir(os.path.join(runAssembly_output, "sff"))
    else:
        for rmopt in [os.path.join(runAssembly_output, opt) for opt in os.listdir(runAssembly_output)]:
            if os.path.isfile(rmopt):
                remove_file(rmopt)
            elif os.path.isdir(rmopt):
                remove_dir(rmopt)
        log.Warning("Data assembly error: Please change the -fc or -sd option.")

    assembly_output = runAssembly_output

    return assembly_output


