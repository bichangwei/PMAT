#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
https://github.com/bichangwei/PMAT

This is part of the PMAT. Based on the results of the assembly 
align with conserved protein-coding genes and select candidate 
seeds. Ensure that blastn is installed and can be called 
directly in the environment variables.
'''

import shutil
import subprocess
import os
import re
import argparse
import time
from log import Log

log = Log()

if os.path.exists(os.path.abspath(os.path.dirname(__file__) + '/../pt_db/Plant_pt.fa')):
    db_path = os.path.join(os.path.abspath(os.path.dirname(__file__) + '/../pt_db'), "Plant_pt.fa")
elif shutil.which("PMAT"):
    db_path = os.path.join(os.path.abspath(os.path.dirname(shutil.which("PMAT")) + '/../pt_db'), "Plant_pt.fa")
else:
    log.Warning("Please check if the pt_db file is installed correctly!")

class SeedFinder:
    def __init__(self, Allcontigs, id_depth, output_path):
        self.Allcontigs = Allcontigs
        self.id_depth = id_depth
        self.output_path = output_path
        self.blastn_out = self._Run_blastn()


    def _Run_blastn(self):
        # Run the blastn.
        # start_time = time.time()

        Blastn_command = ['blastn', '-db', db_path, '-query', self.Allcontigs, '-outfmt', '6', '-num_threads', '30', '-num_alignments', '1', '-max_hsps', '1']
        Blastn_process = subprocess.Popen(Blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Check if child process has terminated. Set and return returncode attribute. Otherwise, returns None.
        # while Blastn_process.poll() is None:
        #     elapsed_time = time.time() - start_time
        #     print(f">>>>>> Find candidate seeds for {elapsed_time:.2f}s <<<<<<", end="\r")
        #     time.sleep(0.1)
        
        # end_time = time.time()
        # elapsed_time = end_time - start_time
        # section_tail(f"\nfinished in {elapsed_time:.2f}s")

        # output the result of blastn and error.
        blastn_out, blastn_err = Blastn_process.communicate(timeout=None) 

        # Error output
        # if blastn_err:
        #     print('\nBLASTn encountered an error:\n' + blastn_err.decode())
        return blastn_out


    def candidate_seeds(self):
        #Find the target contig
        redundant_seed = set()
        if len(self.blastn_out.decode().splitlines()) > 0:
            with open(f'{self.output_path}/PMAT_pt_blastn.txt', 'w') as blt:
                    for line in self.blastn_out.decode().splitlines():
                        blt.write(line+'\n')
                        lines = line.split()
                        # select contigs
                        if int(lines[3]) > 500 and float(lines[2]) > 0.7:
                            redundant_seed.add(re.sub("contig0*", "", lines[0]))
        else:
            log.Warning("No matching candidate contigs found. Please use the graphBuild program to manually modify the --seeds.")
        
        # Sorting according to the contig depth of contig
        candidate_seed = sorted(list(redundant_seed), key=lambda x: float(self.id_depth[str(x)]), reverse=True)
        # log(f"{len(candidate_seed)} contigs are used as candidate seeds")
        return candidate_seed