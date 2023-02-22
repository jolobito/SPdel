## SPDel: Comparing Species delimitation and statistics for DNA Barcoding datasets

SPdel is a pipeline to integrate different single-gene species delimitation methods. SPdel is designed to calculate statistics and compare MOTUs obtained by different methods (e.g. GMYC, PTP, bPTP, BIN, etc).

### Requirement:

* [Anaconda](https://www.anaconda.com/download/)
* [Git](https://git-scm.com/downloads/)

### Installation:

Open Anaconda prompt or terminal:
```
conda env create -f environment.yml
conda activate spdel 
```
### Overview

The pipeline SPdel can run PTP, BPTP and GMYC analyses locally in your computer using original codes and then parse the results. Also, SPdel can interpretate output files from PTP and bPTP web server (https://species.h-its.org/) or python distribution (https://github.com/zhangjiajie/PTP). For GMYC the pipeline works with the outfile from python version (https://github.com/zhangjiajie/pGMYC), but you can use the R version or the web server (https://species.h-its.org/gmyc/) and include your results using -X option. Please check the example files included. 

Please cite original authors of each species delimitation used.

How to Use:

The sequences name should be separate for "_" (e.g. Genus_species_individual) or use -N option for rename sequences

In this example we will use the LA_nominal.fasta dataset included in data folder.

Calculates genetic distances for nominal species
```
SPdel.py path/LA_nominal.fasta -n
```
Calculate PTP, bPTP and GMYC locally and perfom genetic distances analyses
```
SPdel.py path/LA_nominal.fasta -n -P -T -G -t path/LA_nominal.nwk
```
Calculates genetics distances from MOTUs delimited for BIN using outfiles calculated outside the pipeline
```
SPdel.py path/LA_nominal.fasta -B path/LA_BinList.txt -t path/LA_nominal.nwk
```
Compare analyses previously calculated
```
SPdel.py path/LA_nominal.fasta -C n,p,t,g -t path/LA_nominal.nwk
```
Calculate all analyses and compare the results
```
SPdel.py path/LA_nominal.fasta -P -T -G -B path/LA_BinList.txt -t path/LA_nominal.nwk -C n,p,t,g,b
```

Options:   

    -n           For nominal analysis.
    -distance    Substitution model, k for K2p or p for p-distance (default=k)
    -t           Specify the path of the input newick tree for PTP, bPTP and GMYC analysis.
    -N           Specify the path of the text file including the nominal names for rename the sequences.
    -P           Specify the path of the PTP output file.
    -G           Specify the path of the GMYC output file.
    -T           Specify the path of the bPTP output file.     
    -B           Specify the path of the text file including the BIN names obtained from BOLD.
    -X           Specify the path of the text file including the MOTUs names obtained from any external method.
    -D           For Diagnostic character analysis.
    -C           Specify the type of analisys to be compared, include n for nominal, p for PTP, t for bPTP, b for BIN, and any filename  used in X option for external MOTU lists. 
    -code        Specify the genetic code used to test stop codon, VER or INV.

Options for nominal:

    -gen         Position of the genus name in the sequence name when split by "_" (default=1).
    -sp          Position of the species name in the sequence name when split by "_" (default=2)   
    
Options for bPTP
```
    -n_iter       Number of iteration for bPTP analysis (default=10000)
    -sample      Number of sampling for bPTP analysis (default=100)
    -burnin      Burnin for bPTP analysis (default=0.1)     
```

Options for diagnostic character:

    -n_ind       Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3)

### Develop:

After activating the environment

```
conda activate spdel
```

Install the develop requirements:

```
conda install --file requirements_dev.txt
```

 
