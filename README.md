# PMAT - an efficient assembly toolkit for plant mitochondrial genome
<p align="center"><img src="misc/logo.png" alt="PMAT" width="600"></p>

PMAT is an efficient toolkit for mitochondrial genome assembly. It can correct three-generation sequencing data using [canu](https://github.com/marbl/canu) or [NextDenovo](https://github.com/Nextomics/NextDenovo) and then use [Newbler](https://evomics.org/learning/assembly-and-alignment/newbler/) to assemble it. Finally, PMAT can identify mitochondrial genome sequences and generate mitochondrial genome structures.

## Installation

#### recommended installation method
```shell
    wget https://github.com/bichangwei/PMAT/releases/download/v1.1.0/PMAT-1.1.0.tar.gz
    tar -zxvf PMAT-1.1.0.tar.gz
    cd PMAT/bin
    chmod a+x PMAT
    PMAT --help
```

## Repuirement

###### When inputting HiFi reads, the following software are required:
- [**blastn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  Needs to be installed in the `PATH`.
- [**singularity**](https://github.com/YanshuQu/runAssembly)

###### When inputting CLR or ONT reads, the following software are required:

- [**blastn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)     Needs to be installed in the `PATH`.
- [**canu**](https://github.com/marbl/canu) or [**NextDenovo**](https://github.com/Nextomics/NextDenovo)    Suggested to be installed in the `PATH`.
- [**singularity**](https://github.com/YanshuQu/runAssembly)

## options and usage

Run `PMAT --help` to view the program's usage guide.

```
    usage: PMAT <command> <arguments>
    
     ______     ___           __        ____       _____________ 
    |   __  \  |   \        /   |      / __ \     |_____   _____|
    |  |__)  | | |\ \      / /| |     / /  \ \          | |      
    |   ____/  | | \ \    / / | |    / /____\ \         | |      
    |  |       | |  \ \  / /  | |   / /______\ \        | |      
    |  |       | |   \ \/ /   | |  / /        \ \       | |      
    |__|       |_|    \__/    |_| /_/          \_\      |_|


    PMAT            an efficient assembly toolkit for plant mitochondrial genome
    Version         1.1.0
    Contributors    Bi,C. and Han,F.
    Email           bichwei@njfu.edu.cn, hanfc@caf.ac.cn
    
    For more information about PMAT, see https://github.com/bichangwei/PMAT

    optional arguments:
        -h, --help     show this help message and exit
        -v, --version  show program's version number and exit

    Commands:

        autoMito         an efficient assembly toolkit for plant mitochondrial genome
    
        graphBuild  Structure assembly based on Newbler output
```

## Tips

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

## Author
Changwei Bi bichwei@njfu.edu.cn  
Fuchuan Han hanfc@caf.ac.cn

