# <p name="h1">PMAT</p> 
An efficient assembly toolkit for plant mitochondrial genome
<p align="center"><img src="misc/logo.png" alt="PMAT" width="600"></p>

PMAT is an efficient assembly toolkit for assembling plant mitogenomes using third-generation (HiFi/CLR/ONT) sequencing data. PMAT can also be used to assemble chloroplast genomes or animal mitogenomes. 

- [PMAT](#h1)
- [Installation](#C1)
- [Repuirement](#C2)
- [Options and usage](#C3)
    - [autoMito](#C4)
    - [graphBuild](#C5)
- [Examples](#C6)
  - [Demo1](#C6.1)
  - [Demo2](#C6.2)
  - [Demo3](#C6.3)
- [Resulting files](#C7)
- [Version](#C8)
- [Citing PMAT](#C9)

## <a name="C1">Installation </a>

Install by using git
```sh
git clone https://github.com/bichangwei/PMAT.git
cd PMAT/bin
chmod a+x PMAT
PMAT --help
```
Install by downloading the source codes
```sh
wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.2.0.tar.gz
tar -zxvf PMAT-1.2.0.tar.gz
cd PMAT-1.2.0/bin
chmod a+x PMAT
PMAT --help
```

## <a name="C2">Requirement</a>

- [**BLASTn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  Needs to be installed in `PATH`.
- [**Singularity**](https://github.com/YanshuQu/runAssembly)

- [**Canu**](https://github.com/marbl/canu) or [**NextDenovo**](https://github.com/Nextomics/NextDenovo) is required for CLR or ONT sequencing data, which is suggested to be installed in `PATH`.

## <a name="C3">Options and usage</a>

Run `PMAT --help` to show the program's usage guide.
```
usage: PMAT <command> <arguments>

 ______     ___           __        ____       _____________ 
|   __  \  |   \        /   |      / __ \     |_____   _____|
|  |__)  | | |\ \      / /| |     / /  \ \          | |      
|   ____/  | | \ \    / / | |    / /____\ \         | |      
|  |       | |  \ \  / /  | |   / /______\ \        | |      
|  |       | |   \ \/ /   | |  / /        \ \       | |      
|__|       |_|    \__/    |_| /_/          \_\      |_|      

PMAT            An efficient assembly toolkit for plant mitochondrial genome
Version         1.2.0
Contributors    Bi,C. and Han,F.
Email           bichwei@njfu.edu.cn, hanfc@caf.ac.cn

For more information about PMAT, please see https://github.com/bichangwei/PMAT

optional arguments:
-h, --help     show this help message and exit
-v, --version  show program's version and exit

Commands:

    autoMito    One-step de novo assembly of mitochondrial genome. 
                This command can generate the master assembly graph 
                from raw sequencing data directly.

    graphBuild  If 'autoMito' mode fails to generate the mitogenome 
                assembly graph, you can use this command to manually 
                select seeds for assembly.
```

### <a name="C4">autoMito</a>
    
Run `PMAT autoMito --help` to view the usage guide.

```
Required arguments:
  -i INPUT, --input INPUT
                        input raw sequencing file
  -o OUTPUT, --output OUTPUT
                        output directory
  -st SEQTYPE, --seqtype SEQTYPE
                        sequencing platform(ONT/CLR/HiFi)
  -g GENOMESIZE, --genomesize GENOMESIZE
                        Please enter the genome size of the species, such as 1G, 1000M.

optional arguments:
  -h, --help            show this help message and exit
  -tk TASK, --task TASK
                        all/p1/ Default: all
                        all : De novo assembly including error correction for ONT/CLR data and no error correction for HiFi data
                        p1  : Import error-corrected ONT/CLR data for direct assembly
  -cs CORRECTSOFT, --correctsoft CORRECTSOFT
                        Correcting software using nextDenovo or Canu. Default: NextDenovo
  -cp CANU, --canu CANU
                        Please provide the install path of canu.
  -np NEXTDENOVO, --nextDenovo NEXTDENOVO
                        Please provide the install path of nextDenovo.
  -cfg CORRECTCFG, --correctcfg CORRECTCFG
                        config file for nextdenovo correct
  -fc FACTOR, --factor FACTOR
                        Subset extraction of error-corrected ONT, CLR or HiFi data. Sampling ratio factor in 0-1. Default: 1
  -sd SUBSEED, --subseed SUBSEED
                        Sampling set random number seeds, Default: 6
  -bn BREAKNUM, --breaknum BREAKNUM
                        break long reads (>30k) with this. Default: 20000
  -ml MINOVERLAPLEN, --minoverlaplen MINOVERLAPLEN
                        set minimum overlap length. Default: 40
  -mi MINIDENTITY, --minidentity MINIDENTITY
                        set minimum overlap identification. Default: 90
  -cpu CPU              The number of threads. Default: 8
  -l MINLINK, --minLink MINLINK
                        Filter according to the minimum link depth provided by the user
  -v, --version         show program's version and exit
```

**Notes**:
1. Make sure BLASTn was installed in PATH.
2. `-tk` : There are two options for this parameter: 'all' or 'p1'. For ONT or CLR raw data, 'all' is required to correct reads errors and trim the raw data. For error corrected data of ONT/CLR, you can set 'p1' to skip the correct step. For HiFi data, this parameter can be ignored.
3. `-cs` : For ONT or CLR raw data, users should provide the -cs parameter to select the error correction software, default: Nextdenovo.
4. `-cp` : When using Canu for error correction, users need to use -cp parameter to specify the installation path of Canu. This parameter can be ignored when Canu is added to PATH.
5. `-np` : When using NextDenovo to correct errors, users need to use the -np parameter to specify the installation path of NextDenovo. In addition, you need to use canu to trim the data after NextDenovo error correction, so the -cp parameter is needed to specify the installation path of Canu. When NextDenovo and Canu are added to PATH, this parameter can be ignored.
6. `-cfg` : When using NextDenovo error correction, the user needs to specify a config file with -cfg. The contents of the configuration file are recommended to check [NextDenovo](https://nextdenovo.readthedocs.io/en/latest/QSTART.html). Also, it is recommended to add parameter -b for correction_options in config file.
7. `-fc` : This parameter can be used to randomly select a subset of sequencing data for error correction and assembly. Default: all data.
8. `-ml` : Parameters used for assembly, default setting is 40. Recommended setting: 40~200.
9. `-mi` : Parameters used for assembly, default setting is 90. Recommended setting: 90~98.

### <a name="C5">graphBuild</a>

If PMAT fails to generate the assembly graph in 'autoMito' mode, you can use this command to manually select seeds for assembly.

Run `PMAT graphBuild --help` to view the usage guide.

```
Required arguments:
  -c CONTIGGRAPH, --ContigGraph CONTIGGRAPH
                        PMATContigGraph.txt: a file that can get all connections bewteen contigs.
  -a ALLCONTIGS, --AllContigs ALLCONTIGS
                        PMATAllContigs.fna: a file that can get all information of contigs.
  -o OUTPUT, --output OUTPUT
                        output directory
  -gs GENOMESIZE, --genomesize GENOMESIZE
                        Please enter the genome size of the species, such as 1G, 1000M.
  -rs READSIZE, --readsize READSIZE
                        The read size or file for assembly, such as 5G or assembly_seq.cut20K.fasta.

optional arguments:
  -h, --help            show this help message and exit
  -cpu CPU              The number of threads. Default: 8
  -s SEEDS [SEEDS ...], --seeds SEEDS [SEEDS ...]
                        ContigID for extending. Multiple contigIDs should be separated by space. For example: 1 312 356
  -l MINLINK, --minLink MINLINK
                        Filter according to the minimum link depth provided by the user
  -v, --version         show program's version number and exit
```
**Notes**:
1. Make sure BLASTn was installed in PATH.
2. `-c` : PMATContigGraph.txt generated by autoMito command.
3. `-a` : PMATAllContigs.fna generated by autoMito command.
4. `-gs` : The genome size of the species.
5. `-rs` : The amount of data used by the assembly, or provide assembly_seq.cut20K.fasta generated by the graphBuild command.
6. `-s` : Manually select the seeds for the extension, it is recommended to use more than 3 seeds. Use spaces to split between different seed ids, e.g. 1,312,356.

## <a name="C6">Examples</a>

PMAT runtime for different number of threads.

|species|8 cpus|16 cpus|32 cpus|64 cpus|
|-------|------|-------|-------|-------|
|*Arabidopsis thaliana*|13m25.342s|9m29.853s|8m42.429s|9m57.279s|
|*Juncus effusus*|9m37.173s|1m12.433s|4m49.595s|4m40.036s|
|*Malus domestica*|21m12.306s|12m14.663s|7m58.749s|6m48.915s|

**<a name="C6.1">Demo1</a>**

1. Download Arabidopsis thaliana HiFi data
```sh
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Arabidopsis_thaliana_550Mb.fa.gz
```
2. Run autoMito command
```sh
PMAT autoMito -i Arabidopsis_thaliana_550Mb.fa.gz -o ./test1 -st hifi -g 120m
```
3. Run graphBuild (Used when the autoMito command fails to get gfa automatically)
```sh
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -c ./test1/assembly_result/PMATContigGraph.txt -a ./test1/assembly_result/PMATAllContigs.fna -gs 125m -rs ./test1/subsample/assembly_seq.cut20K.fasta -o ./test1_gfa -s 343 345 905 513 1344
```
4. PMAT runtime for different number of threads

```
8 CPUs: 13m25.342s; 16 CPUs: 9m29.853s; 32 CPUs: 8m42.429s; 64 CPUs: 9m57.279s
```

**<a name="C6.2">Demo2</a>**

1. Download Juncus effusus HiFi data
```sh
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Juncus_effusus_216Mb.fa.gz
```
2. Run autoMito command
```sh
PMAT autoMito -i Juncus_effusus_216Mb.fa.gz -o ./test2 -st hifi -g 225m
```
3. Run graphBuild (Used when the autoMito command fails to get gfa automatically)
```sh
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -c ./test2/assembly_result/PMATContigGraph.txt -a ./test2/assembly_result/PMATAllContigs.fna -gs 225m -rs ./test2/subsample/assembly_seq.cut20K.fasta -o ./test2_gfa -s 1 2 457
```
4. PMAT runtime for different number of threads

```
8 CPUs: 9m37.173s; 16 CPUs: 1m12.433s; 32 CPUs: 4m49.595s; 64 CPUs: 4m40.036s
```

**<a name="C6.3">Demo3</a>**

1. Download Malus domestica HiFi data
```sh
wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/Malus_domestica.540Mb.fasta.gz
```
2. Run autoMito command
```sh
PMAT autoMito -i Malus_domestica.540Mb.fasta.gz -o ./test3 -st hifi -g 703m
```
3. Run graphBuild (Used when the autoMito command fails to get gfa automatically)
```sh
# Based on the PMATContigGraph.txt file, manually select 3 or more contigs that match the depth of mitochondrial genome sequencing
PMAT graphBuild -c ./test3/assembly_result/PMATContigGraph.txt -a ./test3/assembly_result/PMATAllContigs.fna -gs 225m -rs ./test3/subsample/assembly_seq.cut20K.fasta -o ./test3_gfa -s 1 2 15391
```
4. PMAT runtime for different number of threads

```
8 CPUs: 21m12.306s; 16 CPUs: 12m14.663s; 32 CPUs: 7m58.749s; 64 CPUs: 6m48.915s
```

## <a name="C7">Resulting files</a>
- The output files include:
  - `*/subsample/assembly_seq_subset.1.0.fasta`, Data for assembly
  - `*/subsample/assembly_seq.cut20K.fasta`, reads trimmed data
  - `*/assembly_result/PMATAllContigs.fna`, The assembly result contains contig sequences
  - `*/assembly_result/PMATContigGraph.txt`, The assembly result contains contig linking relationships
  - `*/assembly_result/PMAT_raw.gfa`, Initial results of mitochondrial genome assembly for visualization
  - `*/assembly_result/PMAT_master.gfa`, Mitochondrial genome assembly optimization results for visualization

## <a name="C8">Version</a>
PMAT version 1.2.0

## <a name="C9">Citing PMAT</a>
Bi C, Shen F, Han F, Qu Y, et al. PMAT: an efficient plant mitogenome assembly toolkit using ultra-low coverage HiFi sequencing data. Unpublished. </br>
Bi C, Qu Y, Hou J, Wu K, Ye Nand Yin T (2022) Deciphering theMulti-Chromosomal MitochondrialGenome of Populus simonii.Front. Plant Sci. 13:914635.doi:10.3389/fpls.2022.914635.
## Author
Changwei Bi, bichwei@njfu.edu.cn  
Fuchuan Han, hanfc@caf.ac.cn

