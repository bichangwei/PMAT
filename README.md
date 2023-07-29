# <p name="h1">PMAT</p> 
An efficient assembly toolkit for plant mitochondrial genome
<p align="center"><img src="misc/logo.png" alt="PMAT" width="600"></p>

PMAT is an efficient assembly toolkit for assembling plant mitogenomes using third-generation (HiFi/CLR/ONT) sequencing data. PMAT can also be used to assemble chloroplast genomes or animal mitogenomes. 

- [PMAT](#h1)
- [Installation instructions](#C1)
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

```sh
git clone https://github.com/bichangwei/PMAT.git
cd PMAT/bin
chmod a+x PMAT
PMAT --help
```

```sh
wget https://github.com/bichangwei/PMAT/archive/refs/tags/v1.2.0.tar.gz
tar -zxvf PMAT-1.2.0.tar.gz
cd PMAT/bin
chmod a+x PMAT
PMAT --help
```

## <a name="C2">Requirements</a>

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

<<<<<<< HEAD
=======
        autoMito         an efficient assembly toolkit for plant mitochondrial genome
    
        graphBuild  Structure assembly based on Newbler output
>>>>>>> f536ea3ebe8563418e096e3fc2b325e5fbf4ba7b
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

<<<<<<< HEAD
## <a name="C6">Example</a>
Tested using data that is already publicly available.
|Species|Data type|SRA/GSA number|
|-------|---------|--------|
|*Pohlia nutans*|HiFi|[CRR383826][CRR383826-data]|
|*Taxus chinensis*|HiFi|[SRR14756467][SRR14756467-data]|
|*Juncus effusus*|HiFi|[ERR8282830][ERR8282830-data]|
|*Arabidopsis thaliana*|HiFi|[CRR302668][CRR302668-data]|
|*Helianthus annuus*|HiFi|[SRR14782853][SRR14782853-data]|
|*Salix wilsonii*|HiFi|[SRR21570388][SRR21570388-data]|
|*Malus domestica*|HiFi|[ERR6939264][ERR6939264-data]|

[CRR383826-data]: https://ngdc.cncb.ac.cn/gsa/browse/CRA006048/CRR383826
[SRR14756467-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR14756467
[ERR8282830-data]: https://www.ncbi.nlm.nih.gov/sra/?term=ERR8282830
[CRR302668-data]: https://ngdc.cncb.ac.cn/gsa/browse/CRA004538/CRR302668
[SRR14782853-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR14782853
[SRR21570388-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR21570388
[ERR6939264-data]:https://www.ncbi.nlm.nih.gov/sra/?term=ERR6939264

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

## <a name="C7">Resulting files</a>
PMAT will generate a series of folders containing the following results:
- ./test*/subsample/assembly_seq.cut20K.fasta </br>> Data for assembly
- ./test*/subsample/assembly_seq.cut20K.fasta </br>> reads trimmed data
- ./test*/assembly_result/PMATAllContigs.fna </br>> The assembly result contains contig sequences
- ./test*/assembly_result/PMATContigGraph.txt </br>> The assembly result contains contig linking relationships
- ./test*/assembly_result/PMAT_raw.gfa </br>> Initial results of mitochondrial genome assembly for visualization
- ./test*/assembly_result/PMAT_master.gfa </br>> Mitochondrial genome assembly optimization results for visualization

## <a name="C8">Version</a>
PMAT version 1.2.0
=======
- **autoMito** -- De novo analysis of sequencing data (hifi, clr and ont).
    
    run `PMAT autoMito --help` view the usage guide for the autoMito command.

    Example:
    ```
    # If canu or NextDenovo is not installed in the PATH when using clr or ont data, you need to provide the parameter -cp or -np.
    For HiFi data    : PMAT autoMito -i hifi.fastq.gz -o hifi_assembly -st hifi -g 500M
    For ONT raw data : PMAT autoMito -i ont.fastq.gz -o ont_assembly -st ont -cfg ont_correct.cfg -tk autoMito -g 500M
    For CLR raw data : PMAT autoMito -i ont.fastq.gz -o clr_assembly -st clr -cfg clr_correct.cfg -tk all -g 500M
    For ONT corrected data : PMAT autoMito -i corrected.fa -o ont_assembly -st ont -tk p1 -g 500M
    For CLR corrected data : PMAT autoMito -i corrected.fa -o clr_assembly -st clr -tk p1 -g 500M
    ```

    notes:

      1. Make sure blastn is installed in the PATH.
      2. If canu or NextDenovo is not installed in the PATH when using clr or ont data, you need to provide the parameter -cp or -np.
        

- **graphBuild** -- Structure assembly based on Newbler results.
    run `PMAT graphBuild --help` view the usage guide for the graphBuild command.

    Example:
    ```shell
    PMAT graphBuild -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs 4G -o output
    PMAT graphBuild -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs hifi.cut20k.fa -s 1 6 8 -o output
    ```
    notes:
    
        When using the graphBuild command, the blastn is required and needs to be installed in the PATH.
    
## Version
PMAT version 1.1.0 (2023.05.05)
>>>>>>> f536ea3ebe8563418e096e3fc2b325e5fbf4ba7b

## <a name="C9">Citing PMAT</a>
Bi C, Shen F, Han F, Qu Y, et al. PMAT: an efficient plant mitogenome assembly toolkit using ultra-low coverage HiFi sequencing data. Unpublished. </br>
Bi C, Qu Y, Hou J, Wu K, Ye Nand Yin T (2022) Deciphering theMulti-Chromosomal MitochondrialGenome of Populus simonii.Front. Plant Sci. 13:914635.doi:10.3389/fpls.2022.914635.
## Author
Changwei Bi, bichwei@njfu.edu.cn  
Fuchuan Han, hanfc@caf.ac.cn

