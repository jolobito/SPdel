SPDel: Comparing Species delimitation and statistics for DNA Barcoding datasets

We developed SPdel a pipeline to integrate different single-gene species delimitation methods. SPdel is designed to calculate statistics and compare MOTUs obtained by different methods (GMYC, PTP, bPTP, BIN, etc).

Installation:

You need to have installed python 3. You can install Anaconda (Python 3) from https://www.anaconda.com/download/ and select the option to put Anaconda in your Path.

Open Anaconda prompt and install dependences:
``` 
conda update -n base conda
conda install numpy
conda install matplotlib
conda install biopython
pip install plotly #conda install -c plotly plotly
```
The pipeline can interpretate output files from PTP and bPTP web server (https://species.h-its.org/) or python distribution (https://github.com/zhangjiajie/PTP). For GMYC the pipeline works with the outfile from python version (https://github.com/zhangjiajie/pGMYC), but you can use R version or web server (https://species.h-its.org/gmyc/) and include your results using -X option. Please check the example files included. 

How to Use:

The sequences name should be separate for "_" (e.g. Genus_species_individual) or use -N option for rename sequences

   
usage: ./SPdel.py path_to_files/ fasta_file -n -P PTP_File -t tree_file -X MOTUList1.txt,MOTUList2.txt -C p,MOTUList1,MOTUList2

usage: ./SPdel.py path_to_files/ fasta_file -n -distance p -code VER

usage: ./SPdel.py path_to_files/ fasta_file -P PTP_File -G GMYC_File -T bPTP_File -t tree_file -B BIN_file


Options:   

    -n           For nominal analysis.
    -distance    Substitution model, k for K2p or p for p-distance (default=k)
    -t           Specify the name of the input newick tree for PTP, bPTP and GMYC analysis.
    -N           Specify the text file including the nominal names for rename the sequences.
    -P           Specify the PTP output file.
    -G           Specify the GMYC output file.
    -T           Specify the bPTP output file.     
    -B           Specify the text file including the BIN names obtained from BOLD.
    -X           Specify the text file including the MOTUs names obtained from any external method.
    -D           For Diagnostic character analysis.
    -C           Specify the type of analisys to be compared, include n for nominal, p for PTP, t for bPTP, b for BIN, and any filename  used in X option for external MOTU lists. 
    -code        Specify the genetic code used to test stop codon, VER or INV.

Options for nominal:

    -gen         Position of the genus name in the sequence name when split by "_" (default=1).

    -sp          Position of the species name in the sequence name when split by "_" (default=2)

Options for diagnostic character:

    -n_ind       Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3)

