# PMAT - an efficient assembly tool for plant mitochondrial genome

## Installation

#### install by using git
    git clone https://github.com/bichangwei/PMAT.git
    cd PMAT/bin
    chmod a+x PMAT
    PMAT --help

#### install from source
    unzip PMAT-main.zip
    cd PMAT.bin
    chmod a+x PMAT
    PMAT --help

## Repuirement
PMAT requires the following software to be installed:

- [**blastn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [**canu**](https://github.com/marbl/canu)
- [**NextDenovo**](https://github.com/Nextomics/NextDenovo)

## options and usage
Run 'PMAT --help' to view the progran's help documentation 

    usage: PMAT <command> <arguments>
    
     ______     ___           __        ____       _____________ 
    |   __  \  |   \        /   |      / __ \     |_____   _____|
    |  |__)  | | |\ \      / /| |     / /  \ \          | |      
    |   ____/  | | \ \    / / | |    / /____\ \         | |      
    |  |       | |  \ \  / /  | |   / /______\ \        | |      
    |  |       | |   \ \/ /   | |  / /        \ \       | |      
    |__|       |_|    \__/    |_| /_/          \_\      |_|


    PMAT            an efficient assembly tool for plant mitochondrial genome
    Version         1.1.0
    Contributors    Bi,C. and Han,F.
    Email           bichwei@njfu.edu.cn, hanfc@caf.ac.cn
    
    For more information about PMAT, see https://github.com/bichangwei/PMAT

    optional arguments:
        -h, --help     show this help message and exit
        -v, --version  show program's version number and exit

    Commands:

        all         an efficient assembly tool for plant mitochondrial genome
    
        graphBuild  Structure assembly based on Newbler output


## Tips

- all -- De novo analysis of sequencing data (hifi, clr and ont)
    
    run 'PMAT all --help' view the help documentation for the command all

    Example:
    
        For HiFi data    : PMAT all -i hifi.fastq.gz -o hifi_assembly -st hifi -g 500M
        For ONT raw data : PMAT all -i ont.fastq.gz -o ont_assembly -st ont -cfg ont_correct.cfg -tk all -g 500M
        For CLR raw data : PMAT all -i ont.fastq.gz -o clr_assembly -st clr -cfg clr_correct.cfg -tk all -g 500M
        For ONT corrected data : PMAT all -i corrected.fa -o ont_assembly -st ont -tk p1 -g 500M
        For CLR corrected data : PMAT all -i corrected.fa -o clr_assembly -st clr -tk p1 -g 500M

    notes:

        When using the all command, canu and NextDenovo software are required. Please save in the
        environment variable or provide the -cp and -np.
        

- graphBuild -- Structure assembly based on Newbler results
    run 'PMAT graphBuild --help' view the help documentation for the command graphBuild

    Example:

        PMAT gfa -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs 4G -o output
        PMAT gfa -c PMATContigGraph.txt -a PMATAllContigs.fna -gs 500M -rs hifi.cut20k.fa -s 1 6 8 -o output

    notes:
    
        when using the graphBuild command, blastn software are required, and needs to be stored in an
        environment variable.
    
## Version
PMAT version 1.10 (2023.05.05)

## Author
Changwei Bi bichwei@njfu.edu.cn
Fuchuan Han hanfc@caf.ac.cn

