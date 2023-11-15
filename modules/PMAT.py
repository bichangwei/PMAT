# !/usr/bin/env python3

"""
https://github.com/bichangwei/PMAT

The script consists of two subroutines, all and graphBuild, with the parameter 'all' for error correction assembly of 
incoming sequencing data and generation of a mitochondrial genome structure graph, and the parameter 'graphBuild' for 
generating the mitochondrial genome structure directly from the assembly results.
"""

import sys
import re
import os
import shutil
import argparse
import time
import gzip
if os.path.dirname(sys.argv[0]):
    if os.path.exists(os.path.abspath(os.path.dirname(sys.argv[0]) + '/../modules/PMAT.py')):
        sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), '../modules'))
    else:
        message = f'Please check if the software is installed correctly in {sys.argv[0]}'
        print('[' + ' \033[1m\033[35mWARNING\033[0m' + ' ] ' + message, file=sys.stderr, flush=True, end='\n')
        sys.exit()
else:
    if shutil.which('PMAT'):
        if os.path.exists(os.path.abspath(os.path.dirname(shutil.which('PMAT')) + '/../modules/PMAT.py')):
            sys.path.append(os.path.join(os.path.dirname(shutil.which('PMAT')), '../modules'))
        else:
            message = f'Please check if the software is installed correctly in {shutil.which("PMAT")}'
            print('[' + ' \033[1m\033[35mWARNING\033[0m' + ' ] ' + message, file=sys.stderr, flush=True, end='\n')
            sys.exit()

import find_condidate_seeds
import find_condidate_seeds_pt
import seeds_extension
import assembly_graph
from check_file import check_file_format, check_file_path, mkdir_file_path
from version import __version__
from log import Log
import get_subsample
import correct_sequences
import break_long_reads
import runassembly
import fastq2fa
import disentangle_mitogenome_from_graph
import gfa2fa
from progressbar import ProgressBar, ProgressBar, Percentage, Bar
# from tqdm import tqdm


def autoMito(args):
    start_time = time.time()
    log.Info(f'PMAT v{__version__}')
    if os.path.exists(f'{args.output}/assembly_result'):
        log.Warning(f'Error: Destination {args.output}/assembly_result already exists.')

    if args.type == "mt":
        mttype = True
        pttype = False
    elif args.type == "pt":
        pttype = True
        mttype = False
    elif args.type == "all":
        mttype = True
        pttype = True
    else:
        log.Warning('Please enter the correct parameters (--type)')

    if args.seqtype.lower() == 'hifi':
        pass
    else:
        if args.canu:
            CANU_PATH = args.canu
        elif shutil.which('canu'):
            CANU_PATH = 'canu'
        else:
            log.Warning('Please check if canu is installed or provide --canu')

    if args.correctsoft:
        if args.correctsoft.lower() == 'nextdenovo':
            if hasattr(args, 'nextDenovo'):
                NEXTDENOVO_PATH = args.nextDenovo
            elif shutil.which('nextDenovo'):
                NEXTDENOVO_PATH = 'nextDenovo'
            else:
                log.Warning('Please check if nextDenovo is installed or provide --nextDenovo')

    mkdir_file_path(args.output)
    Output = os.path.abspath(args.output)
    # mkdir_file_path(f'{Output}/logfile')
    # logfile = f'{Output}/logfile'

    if str(args.genomesize).lower().endswith('g'):
        genomesize = float(str(args.genomesize).lower().replace('g', '')) * 1000000000
    elif str(args.genomesize).lower().endswith('m'):
        genomesize = float(str(args.genomesize).lower().replace('m', '')) * 1000000
    elif str(args.genomesize).lower().endswith('k'):
        genomesize = float(str(args.genomesize).lower().replace('k', '')) * 1000
    else:
        genomesize = float(args.genomesize)


    if args.seqtype.lower() == 'ont':
        readstype = 'nanopore'
    elif args.seqtype.lower() == 'clr':
        readstype = 'pacbio'


    if args.seqtype.lower() == 'hifi':
        high_quality_seq = fastq2fa.fq2fa(args.input)
    else:
        if args.task == 'all':
            log.Info('Run autoMito processes ...')
            if args.correctsoft.lower() == 'nextdenovo':

                if args.correctcfg:
                    high_quality_seq = correct_sequences.ReadsPreprocess(CANU_PATH, args.cpu, readstype, Output, args.correctcfg).NextDenovo_correct(NEXTDENOVO_PATH)
                else:
                    log.Warning('Missing config file, which was required parament when using nextDenovo to correct errors')
            elif args.correctsoft == 'canu':
                if check_file_format(args.input) != 'fastq' and check_file_format(args.input) != 'fasta':
                    log.Warning('Please input the raw data')
                else:
                    log.get_path('Raw data : ' + args.input)
                    high_quality_seq = correct_sequences.ReadsPreprocess(CANU_PATH, args.cpu, readstype, Output).canu_correct(genomesize, args.input)
            pass
        elif args.task == 'p1':
            log.Info('The input file is the corrected ONT/CLR or HiFi data\nif it is more than one file it needs to be merged into one file')
            if check_file_format(args.input) == 'fasta':
                high_quality_seq = args.input
            elif check_file_format(args.input) == 'fastq':
                if args.seqtype.lower() == 'hifi':
                    high_quality_seq = fastq2fa.fq2fa(args.input)
                else:
                    log.Warning(f'Please input corrected ONT/CLR fasta or HiFi reads')
                pass
        else:
            log.Warning('Please enter the correct parameters (--task)')

    assembly_seq_path = judge_subsample(high_quality_seq, Output)
    assembly_seq_cut_path = break_long_reads.BreakLongReads(assembly_seq_path, args.breaknum)
    if args.mem:
        assembly_output = runassembly.run_Assembly(args.cpu, assembly_seq_cut_path, Output, args.minidentity, args.minoverlaplen, mem = args.mem)
    else:
        assembly_output = runassembly.run_Assembly(args.cpu, assembly_seq_cut_path, Output, args.minidentity, args.minoverlaplen)
    file_data_fna_name = f'{assembly_output}/PMATAllContigs.fna'
    file_data_name = f'{assembly_output}/PMATContigGraph.txt'
    assemblysize = os.path.getsize(assembly_seq_path)
    # file_result_name=Output#File name of the final result store

    with open(file_data_name, 'r') as fo:
        data_list = fo.readlines()
    log.section_header("Reading files ...")
    # Result of get_all_connection_from_datafile
    simple_pairs = get_all_connection_from_datafile(data_list)[0]#{'1': 'contig00001',....., '92825': 'contig92825'}
    all_connections = get_all_connection_from_datafile(data_list)[1]  #  [{"68206 3'": "50193 5'"}, .....,{"92534 5'": "92436 3'"}]
    id_depth = get_all_connection_from_datafile(data_list)[2]  # Shorthand ID and sequencing depth for each contig
    id_length = get_all_connection_from_datafile(data_list)[3] # Shorthand ID and sequencing length for each contig
    id_seq = get_all_connection_from_datafile(data_list)[4] # Shorthand ID and sequnece for each contig
    link_depth = get_all_connection_from_datafile(data_list)[5] # Shorthand link depth

    longest_contig = max(id_length, key=id_length.get)
    # longest_contig, longest_len = list(sorted(id_length.items(), key= lambda item: item[1], reverse=True))[0]
    log.Info(f'Contig number  : {len(simple_pairs)}')
    log.Info(f'Longest Contig : contig{longest_contig} {id_length[str(longest_contig)]}bp')
    log.Info(log.dim('-'*60))

    if getattr(args, 'minLink', None):
        minLink = args.minLink
    else:
        minLink = None
            
    # Determine if the user provides parameters for minPath
    # assemblysize = os.path.getsize(assembly_seq_cut_path)
    nucl_contig_depth = assemblysize / genomesize


    if mttype:
        proleptic_connections = []
        mt_condidate_seeds = process_candidate_seeds(simple_pairs, id_length, find_condidate_seeds.SeedFinder, file_data_fna_name, id_depth, Output)
        mt_dynamic_sampleIDs = mt_condidate_seeds
        mt_initial_seeds, mt_main_seeds, raw_seeds,  initial_connections = extend_seeds(mt_condidate_seeds, all_connections, id_depth, 
                                                                                        id_length, link_depth, mt_dynamic_sampleIDs, 
                                                                                        simple_pairs, proleptic_connections, minLink, 
                                                                                        nucl_contig_depth, file_data_fna_name, id_seq, "mt", 2)
        gfa_result(file_data_fna_name, Output, id_length, 
                    id_depth, simple_pairs, initial_connections, 
                    mt_initial_seeds, mt_main_seeds, raw_seeds,'mt')

    if pttype:
        proleptic_connections = []
        pt_condidate_seeds = process_candidate_seeds(simple_pairs, id_length, find_condidate_seeds_pt.SeedFinder, file_data_fna_name, id_depth, Output)
        pt_dynamic_sampleIDs = pt_condidate_seeds
        pt_initial_seeds, pt_main_seeds, raw_seeds ,initial_connections = extend_seeds(pt_condidate_seeds, all_connections, id_depth, 
                                                                                        id_length, link_depth, pt_dynamic_sampleIDs, 
                                                                                        simple_pairs, proleptic_connections, minLink, 
                                                                                        nucl_contig_depth, file_data_fna_name, id_seq, "pt", 10)
        gfa_result(file_data_fna_name, Output, id_length, 
                    id_depth, simple_pairs, initial_connections, 
                    pt_initial_seeds, pt_main_seeds, raw_seeds,'pt')


    log.Info('Generate gfa task end.')
    log.Info('Task over, bye!\n')


    # end_time = time.time()
    # print('Took %f second' % (end_time - start_time))
def graphBuild(args):
    # start_time = time.time()
    if args.type == "mt":
        mttype = True
        pttype = False
    elif args.type == "pt":
        pttype = True
        mttype = False
    elif args.type == "all":
        mttype = True
        pttype = True
    else:
        log.Warning('Please enter the correct parameters (--type)')

    if str(args.genomesize).lower().endswith('g'):
        genomesize = float(str(args.genomesize).lower().replace('g', '')) * 1000000000
    elif str(args.genomesize).lower().endswith('m'):
        genomesize = float(str(args.genomesize).lower().replace('m', '')) * 1000000
    elif str(args.genomesize).lower().endswith('k'):
        genomesize = float(str(args.genomesize).lower().replace('k', '')) * 1000
    else:
        genomesize = float(args.genomesize)

    if os.path.isfile(args.readsize):
        readsize = os.path.getsize(args.readsize)
    else:
        if str(args.readsize).lower().endswith('g'):
            readsize = float(str(args.readsize).lower().replace('g', '')) * 1000000000
        elif str(args.readsize).lower().endswith('m'):
            readsize = float(str(args.readsize).lower().replace('m', '')) * 1000000
        elif str(args.readsize).lower().endswith('k'):
            readsize = float(str(args.readsize).lower().replace('k', '')) * 1000
        else:
            readsize = float(args.readsize)

    file_data_name=args.ContigGraph#Oiginal data file name
    file_data_fna_name=args.AllContigs
    # file_result_name=args.output#File name of the final result store
    mkdir_file_path(f'{args.output}')
    Output = args.output

    with open(file_data_name, 'r') as fo:
        data_list = fo.readlines()
    log.section_header("Reading files ...")
    # Result of get_all_connection_from_datafile
    simple_pairs = get_all_connection_from_datafile(data_list)[0]#{'1': 'contig00001',....., '92825': 'contig92825'}
    all_connections = get_all_connection_from_datafile(data_list)[1]  #  [{"68206 3'": "50193 5'"}, .....,{"92534 5'": "92436 3'"}]
    id_depth = get_all_connection_from_datafile(data_list)[2]  # Shorthand ID and sequencing depth for each contig
    id_length = get_all_connection_from_datafile(data_list)[3] # Shorthand ID and sequencing length for each contig
    id_seq = get_all_connection_from_datafile(data_list)[4] # Shorthand ID and sequnece for each contig
    link_depth = get_all_connection_from_datafile(data_list)[5] # Shorthand link depth

    longest_contig = max(id_length, key=id_length.get)
    # longest_contig, longest_len = list(sorted(id_length.items(), key= lambda item: item[1], reverse=True))[0]
    log.Info(f'Contig number : {len(simple_pairs)}')
    log.Info(f'Longest Contig : contig{longest_contig} {id_length[str(longest_contig)]}bp')
    log.Info(log.dim('-'*60))

    # Determine if the user provides parameters for minLink
    if getattr(args, 'minLink', None):
        minLink = args.minLink
    else:
        minLink = None

    # Result of update_seed_extend
    nucl_contig_depth = readsize / genomesize
    
    if mttype:
        proleptic_connections = []
        if args.seeds:
            mt_condidate_seeds = args.seeds
            log.Info(f'SeedIDs for extending : {str(set(mt_condidate_seeds))}')
            if len(mt_condidate_seeds) == 1:
                for condidate_seed in mt_condidate_seeds:
                    print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="\n")
            else:
                for n, condidate_seed in enumerate(mt_condidate_seeds):
                    if n == 0:
                        print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="  ")
                        time.sleep(0.1)
                    elif n+1 < len(mt_condidate_seeds):
                        print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="  ")
                        time.sleep(0.1)
                    else:
                        print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="\n")
                        time.sleep(0.1)
                pass
            log.Info(log.dim('-'*60))
        else:
            mt_condidate_seeds = process_candidate_seeds(simple_pairs, id_length, find_condidate_seeds.SeedFinder, file_data_fna_name, id_depth, Output)
            
        mt_dynamic_sampleIDs = mt_condidate_seeds
        mt_initial_seeds, mt_main_seeds, raw_seeds,  initial_connections = extend_seeds(mt_condidate_seeds, all_connections, id_depth, 
                                                                                        id_length, link_depth, mt_dynamic_sampleIDs, 
                                                                                        simple_pairs, proleptic_connections, minLink, 
                                                                                        nucl_contig_depth, file_data_fna_name, id_seq, "mt", 2)
        gfa_result(file_data_fna_name, Output, id_length, 
                    id_depth, simple_pairs, initial_connections, 
                    mt_initial_seeds, mt_main_seeds, raw_seeds,'mt')
        
    if pttype:
        proleptic_connections = []
        pt_condidate_seeds = process_candidate_seeds(simple_pairs, id_length, find_condidate_seeds_pt.SeedFinder, file_data_fna_name, id_depth, Output)

        pt_dynamic_sampleIDs = pt_condidate_seeds
        pt_initial_seeds, pt_main_seeds, raw_seeds ,initial_connections = extend_seeds(pt_condidate_seeds, all_connections, id_depth, 
                                                                                        id_length, link_depth, pt_dynamic_sampleIDs, 
                                                                                        simple_pairs, proleptic_connections, minLink, 
                                                                                        nucl_contig_depth, file_data_fna_name, id_seq, "pt", 10)
        gfa_result(file_data_fna_name, Output, id_length, 
                    id_depth, simple_pairs, initial_connections, 
                    pt_initial_seeds, pt_main_seeds, raw_seeds,'pt')


    log.Info('Generate gfa task end.')
    log.Info('Task over, bye!\n')
    # end_time = time.time()
    # print('Took %f second' % (end_time - start_time))


def process_candidate_seeds(simple_pairs, id_length, seed_finder_class, file_data_fna_name, id_depth, Output):
        log.section_header("Candidate seeds search start ...")
        candidate_seeds = seed_finder_class(file_data_fna_name, id_depth, Output).condidate_seeds()
        log.Info(f"{len(candidate_seeds)} contigs are used as candidate seeds")

        for n, condidate_seed in enumerate(candidate_seeds):
        # simple_pairs[str(condidate_seed)] = 'contig' + condidate_seed
            if n == 0:
                print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="  ")
                time.sleep(0.1)
            elif n+1 < len(candidate_seeds):
                print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="  ")
                time.sleep(0.1)
            else:
                print(f'{log.bold_green(simple_pairs[str(condidate_seed)])}: {id_length[str(condidate_seed)]}bp {id_depth[str(condidate_seed)]}X', file=sys.stdout, flush=True, end="\n")
                time.sleep(0.1)
        log.Info(log.dim('-'*60))

        return candidate_seeds

def extend_seeds(condidate_seeds, all_connections, id_depth, id_length, link_depth, pt_dynamic_sampleIDs, simple_pairs, proleptic_connections, min_link, nucl_contig_depth, file_data_fna_name, id_seq, mtpt, dep):
    log.section_header("Seeds extension start ...")
    initial_connections = seeds_extension.Extend_seeds(all_connections, id_depth, id_length, link_depth, 
                                                        pt_dynamic_sampleIDs, simple_pairs, proleptic_connections, 
                                                        min_link, nucl_contig_depth, mtpt).update_seed_extend()  # The list stores contig connections
    log.section_tail("Seeds extension end ...")
    log.Info(log.dim('-'*60))

    initial_seeds = assembly_graph.AssemblyGraph(initial_connections, file_data_fna_name, 
                                                id_length, id_depth, simple_pairs).Target_seeds('init')
    
    raw_add = []
    for candidate_seed in condidate_seeds:
        if str(candidate_seed) not in initial_seeds and float(id_depth[str(candidate_seed)]) > nucl_contig_depth * dep:
            raw_add.append(candidate_seed)

    main_seeds = assembly_graph.AssemblyGraph(initial_connections, file_data_fna_name, id_length, id_depth, simple_pairs).Target_seeds('main')

    # Supplement the missing contig in PMATAllcontigs.txt
    raw_seeds = initial_seeds
    raw_seeds.extend(raw_add)
    add_seq(raw_seeds, file_data_fna_name, simple_pairs, id_seq, id_length)
    return initial_seeds, main_seeds, raw_seeds, initial_connections

def gfa_result(file_data_fna_name, Output, id_length, id_depth, simple_pairs, initial_connections, initial_seeds, main_seeds, raw_seeds,mtpt):
    # Read fna file and output gfa( need more time than Biopython)
    log.Info('Start generating the gfa file ...')
    fna = open(file_data_fna_name, 'r')
    fna_seq = fna.readlines()
    fna.close()
    contig_dict = {}
    start_time = time.time()
    index = 0
    head_seq = None
    while index < len(fna_seq):
        elapsed_time = time.time() - start_time
        if fna_seq[index].startswith('>'):
            head_seq = re.sub('>contig0*', '', fna_seq[index].split()[0].strip())

            if head_seq in initial_seeds:
                contig_dict[head_seq] = ''
        elif head_seq in initial_seeds:
            contig_dict[head_seq] = contig_dict[head_seq] + str(fna_seq[index].strip())
        index = index+1
        # print(f">>>>>> save gfa for {elapsed_time:.2f}s <<<<<<", end="\r")
    end_time = time.time() - start_time
    log.Info(f"save gfa for {end_time:.2f}s ")

    # Output the initial gfa
    mkdir_file_path(f'{Output}/gfa_result')
    init_output = os.path.join(Output, f'gfa_result/PMAT_{mtpt}_raw.gfa')

    if len(raw_seeds) != 0:
        init_gfa = open(init_output, 'w')
        assembly_graph.AssemblyGraph(initial_connections, file_data_fna_name, 
                                    id_length, id_depth, simple_pairs, contig_dict).save_gfa(init_gfa, raw_seeds, 'init')
        init_gfa.close()
        log.Info(f'{str(len(raw_seeds))} contigs are added to a raw graph')
    else:
        log.Error('There is no raw structure for this seeds extension result.')

    # Output the main gfa   
    mkdir_file_path(f'{Output}/gfa_result')
    main_output = os.path.join(Output, f'gfa_result/PMAT_{mtpt}_master.gfa')

    if len(main_seeds) != 0:
        main_gfa = open(main_output, 'w')
        mt_main_connections = assembly_graph.AssemblyGraph(initial_connections, file_data_fna_name, 
                                    id_length, id_depth, simple_pairs, contig_dict).save_gfa(main_gfa, main_seeds, 'main')
        main_gfa.close()
        log.Info(f'{str(len(main_seeds))} contigs are added to a master graph')
        if mtpt == 'mt':
            try:
                log.section_header("PMAT attempting automated loop resolution ...")
                loop_result = disentangle_mitogenome_from_graph.MasterLoops(id_depth, id_length, mt_main_connections, main_output, Output).getloop()
                if loop_result:
                    for num, loop_connection in enumerate(loop_result):
                        num += 1
                        loop_output = os.path.join(Output, f'gfa_result/PMAT_{mtpt}_loop_{num}.gfa')
                        loop_gfa = open(loop_output, 'w')
                        assembly_graph.AssemblyGraph(initial_connections, file_data_fna_name, 
                                                    id_length, id_depth, simple_pairs, contig_dict).save_loop_gfa(loop_gfa, loop_connection)
                        loop_gfa.close()

                        loop_gfa2fa = os.path.join(Output, f'gfa_result/PMAT_{mtpt}_loop_{num}.fasta')
                        gfa2fa.gfa2fa(loop_output, loop_gfa2fa)

                    log.Info("Automated loop resolution completed. Please perform a manual inspection.")
                else:
                    log.Error("Unable to generate cyclic structure.")
            except:
                log.Info("PMAT unable to complete automatic loop resolution.")
    else:
        log.Error('There is no master structure for this seeds extension result.')


def get_all_connection_from_datafile(data_list):
    '''
    Reading the runAssembly results file;
    obtaining information on all contig lengths, depths and connection relationships;
    '''

    all_connections = []  # Stores information about the 5' end and 3' end of each ID  [{"68206 3'": "50193 5'"}, .....,{"92534 5'": "92436 3'"}]   all_connections
    simple_pairs = {}  # The complete ID of each contig and the corresponding shorthand form the dictionary as follows :{'1': 'contig00001',....., '92825': 'contig92825'}
    id_depth = {}  # Shorthand ID and sequencing depth for each contig
    id_length = {} # Shorthand ID and sequnencing length for each contig
    id_seq = [] # Shorthand the id and sequence for contig 
    link_depth = {} # Shorthand the depth of link {('1', '2'): '33', ('2', '6'): '89'}
    
    for line in data_list:

        #Store information of each contig    
        if re.match('[0-9]', line):
            temp_contig = line.strip().split()
            simple_id = str(re.sub('contig0*', '', temp_contig[1]))
            simple_pairs[simple_id] = temp_contig[1]
            id_depth[simple_id] = temp_contig[3]
            id_length[simple_id] = int(temp_contig[2])
            pass

        #Find and store contig connections
        if re.match('C', line):
            connection = {}
            lid = ""
            lio = line.strip().split()
            contig_and_UTR1 = lio[1] + " " + lio[2]
            contig_and_UTR2 = lio[3] + " " + lio[4]  # Put the corresponding contig and edge number into the list
            connection[contig_and_UTR1] = contig_and_UTR2
            all_connections.append(connection)
            link_depth[(lio[1], lio[3])] = lio[5]

        if re.match('I', line):
            match_seq = {}
            lio = line.strip().split()
            match_seq[lio[1]] = lio[2]
            id_seq.append(match_seq)

    return simple_pairs, all_connections, id_depth, id_length, id_seq, link_depth


def add_seq(seeds, file_data_fna_name, simple_pairs, id_seq, id_length):
    '''
    Add the lacking target contig to PMATAllContig.fna
    '''

    with open(file_data_fna_name, 'r') as rfna:
        fna_contend = []    #Set a list that stores contig in PMATAllcontig.fna generated by runAssembly
        for line in rfna.readlines():
            if re.match('>', line):
                fna_contig = re.sub('>', '', line.split()[0])
                fna_contend.append(fna_contig)

    with open(file_data_fna_name, 'a') as wfna:
        for ctg in seeds:
            if simple_pairs[str(ctg)] not in fna_contend: #find the target contig that is lacking in PMATAllContig.fna 
                for match_seq in id_seq:
                    if int(ctg) == int(list(match_seq.keys())[0]):
                        wfna.write(">{} length={}\n{}\n".format(simple_pairs[str(ctg)], id_length[str(ctg)], match_seq[ctg]))


def judge_subsample(corrected_seq, Output):
    assembly_seq_path = get_subsample.subsample(Output, corrected_seq ,args.factor, args.subseed)

    return assembly_seq_path

def get_pmat():
    pmat_art = (
        r"  ______     ___           __        ____       _____________ " + '\n' +
        r" |   __  \  |   \        /   |      / __ \     |_____   _____|" + '\n' +
        r" |  |__)  | | |\ \      / /| |     / /  \ \          | |      " + '\n' +
        r" |   ____/  | | \ \    / / | |    / /____\ \         | |      " + '\n' +
        r" |  |       | |  \ \  / /  | |   / /______\ \        | |      " + '\n' +
        r" |  |       | |   \ \/ /   | |  / /        \ \       | |      " + '\n' +
        r" |__|       |_|    \__/    |_| /_/          \_\      |_|      " + '\n' 
        )

    return pmat_art


if __name__ == '__main__':
    
    log = Log()

    """
    Parse the command line arguments.
    """
    usage = 'PMAT <command> <arguments>'

    description = f"PMAT            an efficient assembly toolkit for plant mitochondrial genome"
    Version =     f"Version         {__version__}"
    Contributor = f"Contributors    Bi,C. and Han,F."
    Email =       f"Email           bichwei@njfu.edu.cn, hanfc@caf.ac.cn"
    github = f"{log.blue('For more information about PMAT, see https://github.com/bichangwei/PMAT')}"

    parser=argparse.ArgumentParser(description=log.purple(get_pmat()) + '\n' + description + 
                                   '\n' + Version + '\n' + Contributor + '\n' + Email + '\n\n' + github, 
                                   formatter_class=argparse.RawTextHelpFormatter, usage=usage)
    parser.add_argument('-v', '--version', action='version', version='PMAT v' + __version__)

    sub_description ="""
    autoMito    One-step de novo assembly of the mitochondrial genome. 
                This command corrects the raw ONT/CLR data or uses 
                the corrected data or HiFi for assembly directly. 
                Based on the assembly result, automatically select 
                seeds for extension and filter false positives to 
                obtain an assembly map of the mitochondrial genome.\n
    graphBuild  If PMAT fails to generate the mitochondrial genome 
                assembly map in one-step assembly, you can use this 
                command by manually select seeds for assembly."""
    subparsers = parser.add_subparsers(title='Commands', description=sub_description, metavar='')

    #--------------------------------------------------------------------
    autoMito_description = f"PMAT: an efficient assembly toolkit for plant mitochondrial genome\n\n{log.blue('For more information about PMAT, see https://github.com/bichangwei/PMAT')}\n\n"
    
    autoMito_example = """Example:
    # If canu or NextDenovo is not installed in the PATH when using ONT/CLR raw data, you need to provide the parameter -cp or -np.
    PMAT autoMito -i hifi.fastq.gz -o hifi_assembly -st hifi -g 500M"""

    autoMito_sub = subparsers.add_parser('autoMito', description=autoMito_description + 
                                         autoMito_example, formatter_class=argparse.RawTextHelpFormatter)

    optional_group=autoMito_sub._action_groups.pop()
    required_group = autoMito_sub.add_argument_group('Required arguments')
    required_group.add_argument('-i', '--input', required=True,
                                help='input raw sequencing file')
    required_group.add_argument('-o', '--output', required=True,
                                help='output directory')
    required_group.add_argument('-st', '--seqtype', required=True,
                                help='sequencing platform(ONT/CLR/HiFi)')
    required_group.add_argument("-g","--genomesize", type=str, required=True, 
                                help='Please enter the genome size of the species, such as 1G, 1000M.')

    optional_group.add_argument('-tk', '--task', required=False, default='all',
                                help='all/p1/ Default: all\nall : De novo assembly including error correction for ONT/CLR data and no error correction for HiFi data\np1  : Import error-corrected ONT/CLR data for direct assembly')
    optional_group.add_argument('-tp', '--type', required=False, default='mt',
                                help='mt/pt/all Default: mt\nmt   : Assembling the mitochondrial genome\npt   : Assembling the chloroplast genome\nall  : Assembling the mitochondrial and chloroplast genomes')
    optional_group.add_argument('-cs', '--correctsoft', required=False, default='nextDenovo',
                                help='Correcting software using nextDenovo or Canu. Default: nextDenovo')
    optional_group.add_argument('-cp', '--canu', required=False,
                                help='Please provide the install path of canu.')
    optional_group.add_argument('-np', '--nextDenovo', required=False,
                                help='Please provide the install path of nextDenovo.')
    optional_group.add_argument('-cfg', '--correctcfg', required=False,
                                help='config file for nextdenovo correct')
    optional_group.add_argument('-fc', '--factor', type=float, required=False, default=1,
                                help='Subset extraction of error-corrected ONT, CLR or HiFi data. Sampling ratio factor in 0-1. Default=1')
    # optional_group.add_argument('-rn', '--random', type=int, required=False, default=6,
    #                             help='Random number seeding when extracting subsets. Default:6')
    optional_group.add_argument('-sd', '--subseed', required=False, default=6,
                                help='Sampling set random number seeds, Default=6')
    # optional_group.add_argument('-tl', '--trimlen', type=int, required=False, default=500,
    #                             help='ignore reads with length < this. Default: 500')
    optional_group.add_argument('-bn', '--breaknum', type=int, required=False, default=20000,
                                help='break long reads (>30k) with this. Default: 20000')
    optional_group.add_argument('-ml', '--minoverlaplen', type=int, required=False, default=40,
                                help='set minimum overlap length. Default:40')
    optional_group.add_argument('-mi', '--minidentity', type=int, required=False, default=90,
                                help='set minimum overlap identification. Default:90')
    optional_group.add_argument('-cpu',type=int, required=False, default=8,
                                help='Setting the number of threads. Default: 8')
    # optional_group.add_argument("-s","--seeds",required=False, nargs='+',
    #                             help='SeedIDs for extending. Multiple parameters should be separated by Spaces. For example: 1 312 356')
    # optional_group.add_argument("-l","--minLen", type=int, required=False, default=100, 
    #                             help='Filter according to the minimum contig length provided by the user')
    # optional_group.add_argument("-p","--minPath", type=int, required=False, 
    #                             help='Filter according to the minimum path depth provided by the user')
    optional_group.add_argument("-l","--minLink", type=float, required=False, 
                                help='Filter according to the minimum link depth provided by the user')
    optional_group.add_argument("-m", "--mem", action="store_true", required=False,
                                help='Flag to keep sequence data in memory to speed up cup time')
    optional_group.add_argument('-v', '--version', action='version', version='PMAT v' + __version__,)
    autoMito_sub._action_groups.append(optional_group)
    autoMito_sub.set_defaults(func=autoMito)

    #--------------------------------------------------------------------------
    graphBuild_description = f"Structure assembly based on runAssembly output\n\n{log.blue('For more information about PMAT, see https://github.com/bichangwei/PMAT')}\n\n"
    graphBuild_example = """Example:
    PMAT graphBuild -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs 4G -o output
    PMAT graphBuild -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs assembly_seq.cut20K.fasta -s 1 6 8 -o output"""
    
    graphBuild_sub = subparsers.add_parser('graphBuild', description=graphBuild_description + graphBuild_example, formatter_class=argparse.RawTextHelpFormatter)
    
    optional_group=graphBuild_sub._action_groups.pop()

    required_group = graphBuild_sub.add_argument_group('Required arguments')
    required_group.add_argument("-c","--ContigGraph", type=str, required=True, 
                                help='PMATContigGraph.txt: a file that can get all connections between contigs.')
    required_group.add_argument("-a","--AllContigs", type=str, required=True, 
                                help='PMATAllContigs.fna: a file that can get all the information of contigs.')
    required_group.add_argument('-o', '--output', required=True,
                                help='output directory')
    required_group.add_argument("-gs","--genomesize", type=str, required=True, 
                                help='Please enter the genome size of the species, such as 1G, 1000M.')
    required_group.add_argument("-rs", "--readsize",type=str, required=True,
                                help='The read size or file for assembly, such as 5G or assembly_seq.cut20K.fasta.')

    optional_group.add_argument('-tp', '--type', required=False, default='mt',
                                help='mt/pt/all Default: mt\nmt :Assembling the mitochondrial genome\npt  : Assembling the chloroplast genome\nall  :Assembling the mitochondrial and chloroplast genomes')
    optional_group.add_argument('-cpu',type=int, required=False, default=8,
                                help='cpus to use. Default: 8')
    optional_group.add_argument("-s","--seeds",required=False, nargs='+',
                                help='ContigID for extending. Multiple contigIDs should be separated by space. For example: 1 312 356')
    # optional_group.add_argument("-l","--minLen", type=int, required=False, default=100, 
    #                             help='Filter according to the minimum contig length provided by the user')
    #***# optional_group.add_argument("-p","--minPath", type=int, required=False, 
    #                             help='Filter according to the minimum path depth provided by the user')
    optional_group.add_argument("-l","--minLink", type=float, required=False, 
                                help='Filter according to the minimum link depth provided by the user')
    optional_group.add_argument('-v', '--version', action='version', version='PMAT v' + __version__,)
    graphBuild_sub._action_groups.append(optional_group)
    graphBuild_sub.set_defaults(func=graphBuild)

    #----------------------------------------------------------------------------------------

    args=parser.parse_args()

    if len(vars(args)) == 0:
        parser.print_help()
        sys.exit(1)
    else:
        args.func(args)
        if hasattr(args, 'autoMito'):
            autoMito(args)
        elif hasattr(args, 'graphBuild'):
            graphBuild(args)
