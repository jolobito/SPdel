SPDel: Comparing Species delimitation and statistics for DNA Barcoding datasets

We developed SPdel a pipeline to integrate different single-gene species delimitation methods. SPdel is designed to calculate and compare MOTUs obtained by the different methods (GMYC, PTP, bPTP, Spider).

Installation:

You need to have installed python 3. You can install Anaconda (Python 3) from https://www.anaconda.com/download/#windows and select the option to put Anaconda in your Path.

Open Anaconda prompt and install dependences:
 
conda update -n base conda
conda install setuptools numpy
conda install scipy lxml matplotlib
conda install mysql-python # conda install -c anaconda mysql-python
conda install biopython
pip install plotly # conda install -c plotly plotly

How to Use:

The sequences name should be separate for "_" (e.g. Genus_species_individual) or use -N option for rename sequences

usage: ./SPdel.py path_to_files/ fasta_file -a nptgs -T tree_file -X SC_MOTUList1.txt,SC_MOTUList2.txt -C n,p,t,g,s,MOTUList1,SC_MOTUList2 -code VER
usage: ./SPdel.py path_to_files/ fasta_file -a n -distance p
usage: ./SPdel.py path_to_files/ fasta_file -a nptgbs -T tree_file -B BIN_file
usage: ./SPdel.py path_to_files/ fasta_file -a ptg -C p,t,g

Options:
    -a           Specify the type of analisys to be performed, include n for nominal, p for PTP, t for bPTP, b for BIN, s for Spider, c for Compare, x for external MOTU list, and d for Diagnostic character.

    -distance    Substitution model, k for K2p or p for p-distance (default=k)
    -T           Specify the name of the input newick tree for PTP, bPTP and GMYC analysis.

    -N           Specify the text file including the nominal names for rename the sequences.

    -B           Specify the text file including the BIN names obtained from BOLD, mandatory b option.

    -X           Specify the text file including the MOTUs names obtained from any external method, mandatory x option.

    -code        Specify the genetic code used to test stop codon, VER or INV.

Options for nominal:

    -gen         Position of the genus name in the sequence name when split by "_" (default=1).

    -sp          Position of the species name in the sequence name when split by "_" (default=2)

Options for bPTP:

    -n_iter       Number of iteration for bPTP analysis (default=10000)

    -sample      Number of sampling for bPTP analysis (default=100)

    -burnin      Burnin for bPTP analysis (default=0.1)

Options for diagnostic character:

    -n_ind       Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3)

