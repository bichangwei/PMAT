#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
https://github.com/bichangwei/PMAT

Identification of the three generations of input data (ONT, CLR, HiFi), 
if the input data is HiFi then directly assembled, otherwise need to 
correct errors first.
'''

import configparser
import subprocess
import os 
from log import Log
from check_file import mkdir_file_path

log = Log()

class ReadsPreprocess:
    log = Log()
    def __init__(self, canu_path=None, cpu=None, readstype=None, output_path=None, cfg=None):
        # self.logfile = logfile
        self.canu_path = canu_path
        self.cpu = cpu
        self.readstype = readstype
        self.output_path = output_path
        self.cfg = cfg



    def canu_correct(self, genomeSize, seq_path):
        log.section_header("Reads correct start...")
        mkdir_file_path(f'{self.output_path}/correct_out')
        command = f'{self.canu_path} -correct -p PMAT -d {self.output_path}/correct_out genomeSize={genomeSize} useGrid=false batThreads={self.cpu} -{self.readstype} {seq_path}'.split(' ')
        
        # canulog = open(f'{self.logfile}/canu_assembly.log', 'w')
        canu_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            if canu_process.poll() is not None:
                break
            else:
                # canulog.write(canu_process.stdout.readline().decode().strip() + '\n')
                log.log(canu_process.stdout.readline().decode().strip())
        canu_process.communicate()
        # canulog.close()
        
        log.section_tail("Reads correct end.")
        log.get_path(f'Reads correct result : {self.output_path}/correct_out')

        canu_corrected_seq = f'{self.output_path}/correct_out/PMAT.correctedReads.fasta.gz'
        self.canu_trim(genomeSize, canu_corrected_seq)

        high_quality_seq = f'{self.output_path}/trim_out/PMAT.trimmedReads.fasta.gz'
        return high_quality_seq

    def NextDenovo_correct(self, nextDenovo_path):
        '''
        nextDenovo run.cfg
        '''

        if self.config_info()[1] == 'correct':
            log.section_header("Reads correct start...")
        else:
            log.section_header("Start using nextDenovo to correct and assemble...\n")

        command = [nextDenovo_path, self.cfg]
        # mkdri_file_path(f'{self.output_path}/correct_out/')
        workdir = self.config_info()[0]
        # nextdenovolog = open(f'{self.logfile}/nextdenovo_correct.log', 'w')
        nextdenovo_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:
            if nextdenovo_process.poll() is not None:
                break
            else:
                # nextdenovolog.write(nextdenovo_process.stdout.readline().decode().strip() + '\n')
                log.log(nextdenovo_process.stdout.readline().decode().strip())
        # Wait for the process to finish
        nextdenovo_process.communicate()
        # nextdenovolog.close()
        

        if self.config_info()[1] == 'correct':
            log.section_tail("Reads correct end...")
        else:
            log.section_tail("nextDenovo for correction and assembly end...\n")

        cns_path = f'{self.config_info()[0]}/02.cns_align/01.seed_cns.sh.work/'

        with open(f'{self.output_path}/PMAT_correctedReads.fasta', 'w') as fcns:
            # Obtain the result path of nextDenovo error correction
            for root, dirs, files in os.walk(cns_path):
                for file in files:
                    if file == 'cns.fasta':
                        # Obtain the absolute path of cns.fasta
                        cns_file = os.path.join(root, file)
                        with open(cns_file) as fin:
                            # merge all cns.fasta
                            fcns.write(fin.read())

        nextdenovo_corrected_seq = f'{self.output_path}/PMAT_correctedReads.fasta'
        self.canu_trim(self.config_info()[2], nextdenovo_corrected_seq)

        high_quality_seq = f'{self.output_path}/trim_out/PMAT.trimmedReads.fasta.gz'
        return high_quality_seq


    def canu_trim(self, genomeSize, seq_path):
        '''
        Trim the output of the correction
        '''

        log.section_header("Reads trim start...")
        mkdir_file_path(f'{self.output_path}/trim_out')
        command = f'{self.canu_path} -trim -p PMAT -d {self.output_path}/trim_out genomeSize={genomeSize} useGrid=false batThreads={self.cpu} -corrected -{self.readstype} {seq_path}'.split(' ')
        
        # trimlog = open(f'{self.logfile}/canu_trim.log', 'w')
        trim_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:
            if trim_process.poll() is not None:
                break
            else:
                log.log(trim_process.stdout.readline().decode().strip())

        trim_process.communicate()
        # trimlog.close()

        log.section_tail("Reads trim end.")
        log.get_path(f'Reads trim result : {self.output_path}/trim_out')

    
    def config_info(self):
        run_cfg = configparser.ConfigParser()
        run_cfg.read(self.cfg)
        workdir = run_cfg['General']['workdir']
        abs_workdir = os.path.join(os.path.dirname(os.path.abspath(self.cfg)), workdir)
        task = run_cfg['General']['task']
        genomeSize = run_cfg['correct_option']['genome_size']
        return abs_workdir, task, genomeSize
    
    
if __name__ == '__main__':
    # ReadsPreprocess('/home/hanfc/software/canu-2.2/canu-2.2/bin/canu', '30', 'pacbio', './canu_test').canu_correct('100M', '/home/hanfc/pub/Linux_002/workspace/Lour_data/PicBio/sample/SCZ/m54059_161226_153304.subreads.fq.gz /home/hanfc/pub/Linux_002/workspace/Lour_data/PicBio/sample/m54059_161221_221131.subreads.fq.gz')
    ReadsPreprocess('/home/hanfc/software/canu-2.2/canu-2.2/bin/canu', '30', 'nanopore', '/home/hanfc/pub/Linux_001/workspace/mitogenome/Litsea_cubeba/Canu_newbler/script/debbug/nextDenovo_test', './test_data/run.cfg').NextDenovo_correct('/home/hanfc/software/NextDenovo/nextDenovo')
    # ReadsPreprocess('/home/hanfc/software/canu-2.2/canu-2.2/bin/canu', '30', 'pacbio', './canu_test').canu_trim('1.2G', '/home/hanfc/pub/Linux_001/workspace/mitogenome/Litsea_cubeba/Canu_newbler/1.canu_output/Lour.correctedReads.10G.fasta')