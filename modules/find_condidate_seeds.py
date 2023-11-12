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

if os.path.exists(os.path.abspath(os.path.dirname(__file__) + '/../Conserved_PCGs_db/Plant_conserved_mtgene_nt.fa')):
    db_path = os.path.join(os.path.abspath(os.path.dirname(__file__) + '/../Conserved_PCGs_db'), "Plant_conserved_mtgene_nt.fa")
elif shutil.which("PMAT"):
    db_path = os.path.join(os.path.abspath(os.path.dirname(shutil.which("PMAT")) + '/../Conserved_PCGs_db'), "Plant_conserved_mtgene_nt.fa")
else:
    log.Warning("Please check if the Conserved_PCGs_db file is installed correctly!")

class SeedFinder:
    def __init__(self, Allcontigs, id_depth, output_path):
        self.Allcontigs = Allcontigs
        self.id_depth = id_depth
        self.output_path = output_path
        self.conserved_PCGs = "atp1 atp4 atp6 atp8 atp9 cob cox1 cox2 cox3 mttB matR rpl2 rpl10 rpl16 rps1 rps3 rps4 rps7 rps12 nad1 nad2 nad3 nad4 nad4L nad5 nad6 nad7 nad9 ccmB ccmC ccmFn ccmFc sdh3 sdh4".split(' ')
        self.PCGs_exon = "rpl2 rps3 nad1 nad2 nad4 nad5 nad7 ccmFc".split(' ')
        self.PCGs_len = self._PCGs_length()
        self.blastn_out = self._Run_blastn()


    def _PCGs_length(self):
        fna = open(db_path, 'r')
        PCGs_line = fna.readlines()
        fna.close()
        
        PCGs_seq = {}
        PCGs_len = {}
        index = 0
        # print(PCGs_line)
        while index < len(PCGs_line):
            if re.match('>', PCGs_line[index]):
                PCGs_seq[PCGs_line[index]] = PCGs_line[index+1]
            index = index+1

        for PCGs in self.conserved_PCGs:
            list_length = []

            for PCG_id, PCG_seq in PCGs_seq.items():
                if PCGs == str(PCG_id.strip().split('_')[2]):
                    list_length.append(len(PCG_seq))

            PCGs_len[PCGs] = int(sum(list_length) / len(list_length))

        return PCGs_len


    def _Run_blastn(self):
        # Run the blastn.
        # start_time = time.time()

        Blastn_command = ['blastn', '-db', db_path, '-query', self.Allcontigs, '-outfmt', '6', '-num_threads', '30', '-num_alignments', '1', '-max_hsps', '1']
        Blastn_process = subprocess.Popen(Blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Check if child process has terminated. Set and return returncode attribute. Otherwise, returns None.
        # while Blastn_process.poll() is None:
        #     elapsed_time = time.time() - start_time
        #     print(f">>>>>> Find condidate seeds for {elapsed_time:.2f}s <<<<<<", end="\r")
        #     time.sleep(0.1)
        
        # end_time = time.time()
        # elapsed_time = end_time - start_time
        # section_tail(f"\nfinished in {elapsed_time:.2f}s")

        # output the result of blastn and error.
        blastn_out, blastn_err = Blastn_process.communicate(timeout=None) 

        # Error output
        if blastn_err:
            print('\nBLASTn encountered an error:\n' + blastn_err.decode())
        return blastn_out


    def condidate_seeds(self):
        #Find the target contig
        PCGs_len = self.PCGs_len
        blast_info = [] # [['atp1', 'contig00001', '90.9', '501'], ['apt2', 'contig00002', '891', '1030'] ...]
        redundant_seed = set()
        with open(f'{self.output_path}/PMAT_mt_blastn.txt', 'w') as blt:
            for line in self.blastn_out.decode().splitlines():
                blt.write(line+'\n')
                lines = line.split()
                PCGs = re.sub(r".*_", "", lines[1])
                PCGs = re.sub(r"-.*", "", PCGs)
                # select contigs
                if PCGs in self.PCGs_exon and int(lines[3]) > 500 and float(lines[2]) > 0.85:
                    redundant_seed.add(re.sub("contig0*", "", lines[0]))
                    
                elif PCGs in PCGs_len.keys() and int(lines[3]) > float(PCGs_len[PCGs]) * 0.9 and float(lines[2]) > 0.95:
                    blast_info.append([PCGs, lines[0], lines[2], lines[3]])

        PCG_contigs = {} # {'apt1': ("contig00001",“”contig00007")， 'atp2': ("contig000012", "contig000003")...}
        blast_ident = {} # {"contig00001":"90.9", "contig00002":"89.1", "contig000010":"99.3"...}
        blast_len = {} # {"contig00001":"501", "contig00002":"1030", "contig000010":"2010"...}
        contig_PCGs = {} # {"contig00001":"apt1", "contig00002":"atp1", "contig000010":"nad9" ...}

        for PCGs in PCGs_len.keys():
            PCGs_map = set()
            for blast_match in blast_info:
                if PCGs == blast_match[0]:
                    PCGs_map.add(blast_match[1])
                    contig_PCGs[blast_match[1]] = blast_match[0]
                    blast_ident[blast_match[1]] = blast_match[2]
                    blast_len[blast_match[1]] = blast_match[3]
            if bool(PCGs_map):
                PCG_contigs[PCGs] = PCGs_map
        
        for PCGs, contigs in PCG_contigs.items():
            if len(contigs) > 1:
                max_len = max(contigs, key=lambda x: blast_len[x])
                max_ident = max(contigs, key=lambda x:blast_ident[x])
                if max_len == max_ident:
                    redundant_seed.add(re.sub(r"contig0*", "", max_len))
                else:
                    redundant_seed.add(re.sub(r"contig0*", "", max_len))
                    redundant_seed.add(re.sub(r"contig0*", "", max_ident))
            elif len(contigs) == 1:
                redundant_seed.add(re.sub(r"contig0*", "", contigs.pop()))
        # Sorting according to the contig depth of contig
        condidate_seed = sorted(list(redundant_seed), key=lambda x: float(self.id_depth[str(x)]), reverse=True)
        # log(f"{len(condidate_seed)} contigs are used as candidate seeds")
        return condidate_seed

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description="Used to view blastn output results")
    parser.add_argument("--AllContigs","-a",help='454AllContigs.fna: a file that can get all information of contigs.',type=str,required=True)
    parser.add_argument("--Output", "-o", help='Output file for Blastn', type=str,required=False,default="blastn_conversed_PCGs.out")
    args=parser.parse_args()

    Allcontigs=args.AllContigs
    id_depth = {}

    with open(Allcontigs, 'r') as fo:
        for line in fo.readlines():
            #Store information of each contig
            if re.match('[0-9]', line):
                temp_contig = line.strip().split()
                simple_id = str(re.sub('contig0*', '', temp_contig[1]))
                id_depth[simple_id] = temp_contig[3]

    fi = open(args.Output, 'w')
    fi.write("# Fields: query\tsubject\tidentity\tlength\tmismatchs\tgaps\tq.start\tq.end\ts.start\ts.end\tevalue\tscore\n")
    for line in SeedFinder(Allcontigs, id_depth).blastn_out.decode().splitlines():
        fi.write(line+'\n')
    fi.close()
