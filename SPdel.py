#! /usr/bin/env python

# SPdel v 2.0
import os
import sys
import logging
# from sys import platform as sys_pf

try:
    # if sys_pf == 'darwin':
    #     import matplotlib

    #     matplotlib.use("Qt5Agg")

    # from matplotlib import pyplot as plt
    from Bio import SeqIO, Phylo, Seq
    from spdelib import Diagnoser
    from spdelib import Matrian
    import toytree
    import toyplot
    import toyplot.svg
    import toyplot.pdf
    import pandas as pd

except ImportError:
    sys.stderr.write("Please install all dependences first.\n")
    sys.stderr.flush()
    sys.exit()

def check_MOTUlist(fasta, MOTUList):
    """A function to check if the MOTU list and fasta file have same individuals.
    Parameters
    ----------
    fasta: str
        The fasta file.
    fileList: str
        The MOTU list file.        
    Returns
    -------
    True if the files have same numeber of individuals.
    """
    # with open(fasta, newline='') as fasta:
    fasta.seek(0)
    handle = fasta
    seqnames=[]
    n_seq = 0
    for seq in SeqIO.parse(handle, "fasta"):
        n_seq += 1
        seqnames.append(seq.name)                            
    namesfile = MOTUList.index.tolist()
    names= set.intersection(set(seqnames), set(namesfile))    
    if len(names) != len(seqnames):        
        logging.info('Error: Different number of individuals en MOTU list and fasta')
        logging.info('Names found only in fasta:'+ ', '.join(set(seqnames)-set(namesfile))) #extra
        logging.info('Names found only in MOTUList:'+ ', '.join(set(namesfile)-set(seqnames))) #extra
        return False
 
class sorting_data:
    """Initial check for fasta
    
    Parameters
    ----------
    basepath: str
        The folder for outputs.    
    fasta: str
        The fasta file.
    typeCODE: str
        Type of genetic code used to check stop codon, use 'VER' or 'INV'      
    outname: str
        The name of fasta sorted by name.   
    """
    def __init__(self, basepath, fasta, typeCODE=None, outname='Nominal/sorted.fasta',tree=None):
        self.basepath = basepath
        self.fasta = open(fasta, newline='')
        self.outname = outname
        self.typeCODE = typeCODE
        self.tree = tree

    def check_fasta(self):
        """To check if sequences are aligned."""
        self.fasta.seek(0)
        handle = self.fasta
        n = len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
        for seq in SeqIO.parse(handle, "fasta"):
            if (len(seq)) == n:
                continue
            else:
                return False
        return True

    def describe_fasta(self):
        """To print the main characteristics of fasta file."""        
        self.fasta.seek(0)
        handle = self.fasta
        n = len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
        n_seq = 0
        for seq in SeqIO.parse(handle, "fasta"):
            n_seq += 1
        logging.info('Fasta file with ' + str(n_seq) + ' sequences and ' + str(n) + ' base pairs')

    def stop_codon(self): # improve to return the position in the ORF? add N to made multiple of three
        """To check if there are at least one ORF without stop codon"""    
        logging.info('Checking stop codons using ' + str(self.typeCODE) + ' genetic code')
        if self.typeCODE == 'VER':
            table_code=2
        elif self.typeCODE == 'INV':
            table_code=5
        else:
            table_code = self.typeCODE
        allcheck = True
        self.fasta.seek(0)
        handle = self.fasta
        ORF = [0, 1, 2]
        for seq in SeqIO.parse(handle, "fasta"):
            seq.seq=Seq.Seq(str(seq.seq).replace('-',''))
            seq_rc=seq.reverse_complement()
            codonPass = []
            for i in ORF:
                trans=seq[i:].translate(table=table_code)
                trans_rc=seq_rc[i:].translate(table=table_code)
                if '*' in trans:
                    codonPass.append(1)
                else:
                    codonPass.append(0)
                if '*' in trans_rc:
                    codonPass.append(1)
                else:
                    codonPass.append(0)
            # logging.info(codonPass)
            if 0 in codonPass:
                pass
            else:
                allcheck = False
                logging.info('WARNING: No ORF without stop codon were found in sequence ' + seq.id)
        if allcheck == True:
            logging.info('All sequences have at least one ORF without stop codons')
            
    def sort_fasta(self):
        """Sort the fasta file by name in sequences."""    
        self.fasta.seek(0)
        tmp = os.path.join(self.basepath,self.outname)
        file_sorted = open(tmp, 'w+')        
        handle = self.fasta
        fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        records=[]
        for k in sorted(fasta_dict):
            records.append(fasta_dict[k])
        SeqIO.write(records,file_sorted,'fasta')        
        file_sorted.close()
        return file_sorted

    def check_pattern(self):
        """Check if the patter genus_species_ID in the sequences name."""    
        self.fasta.seek(0)
        for line in self.fasta:
            if line.startswith('>'):
                pattern = line.split('_')
                if len(pattern) < 3:
                    return False

    def check_tree(self, treetype='newick'):
        """Check if the names in the tree are exactly the same than the fasta
        Parameters
        ----------
        fasta: str
            The sorted fasta file.
        tree: str
            The tree file.        
        """
        fasta = SeqIO.parse(self.fasta,'fasta')
        tree = Phylo.read(self.tree,treetype)
        seqnames=[]
        treenames=[]
        for seq in fasta:
            seqnames.append(seq.name)
        for leaf in tree.get_terminals():
            treenames.append(leaf.name)
        names= set.intersection(set(seqnames), set(treenames))    
        if len(names) != len(seqnames):
            logging.info('Different number of individuals in tree and fasta')
            logging.info('Names found only in fasta:'+ ', '.join(set(seqnames)-set(treenames))) #extra
            logging.info('Names found only in tree:'+ ', '.join(set(treenames)-set(seqnames))) #extra
            return False

def MOTU_listPTP(path,infile,analisis):
    """Read the input precalclate and construct a df with the MOTUs
    Parameters
    ----------
    infile: str
        species delimitation results file precalculated.
    Returns
    -------
    A df with MOTUs
    """          
    listPTP = open(os.path.join(path, infile))
    listPTP.seek(0)
    name = []
    ind = []
    n=0
    for line in listPTP:
        line = line.strip()
        if not line.strip(): continue
        if line.startswith("#"):
            pass
        elif line.startswith("Species"):
            n+=1
        else:
            for i in line.split(','):
                name.append('MOTU_'+str(n))
                ind.append(i)
    datas = {analisis: name}
    PTP_df = pd.DataFrame(datas, index=ind)  
    return PTP_df


def ptp(path,tree):
    """Run the PTP analysis"""
    from spdelib import PTP     
    PTP.main(['-t', tree, '-o', os.path.join(path, 'PTP', 'PTP'), '-r'])
    
def bptp(path,tree,niter,sample,burnin):
    """Run the bPTP analysis"""
    from spdelib import bPTP          
    bPTP.main(['-t', tree, '-o', os.path.join(path, 'bPTP', 'bPTP'), '-s', '1234', '-r', '-i', niter, '-n',
                sample, '-b', burnin])

def gmyc_run(path,tree):
    """Run the GMYC analysis""" 
    from spdelib import GMYC
    GMYC.main(['-t', tree, '-ps','path',path])

def MOTU_listGMYC(path,infile,analisis):
    """Read the input precalclate and construct a df with the MOTUs
    Parameters
    ----------
    infile: str
        species delimitation results file precalculated.
    Returns
    -------
    A df with MOTUs
    """          
    listGMYC = open(os.path.join(path, infile))
    listGMYC.seek(0)
    name = []
    ind = []
    n=0
    for line in listGMYC:
        line = line.strip()
        line = line.replace(" ", "")
        line = line[:-1]
        if not line.strip(): continue
        elif line.startswith("Species"):
            n+=1
        else:
            for i in line.split(','):
                name.append('MOTU_'+str(n))
                ind.append(i)
    datas = {analisis: name}
    GMYC_df = pd.DataFrame(datas, index=ind)  
    return GMYC_df


def MOTU_X(path, fasta, folder, MOTUDict, k):
    """Rename the fasta name sequences using MOTUs Dictionary
    Parameters
    ----------
    path: str
        The path to folder with the input file.
    fasta: str
        The name of sorted fasta file.
    folder: str
        The folder name for outputs
    MOTUDict: dict
        Dictionary with species delimitation results file precalculated.        
    """    
    if not os.path.exists(os.path.join(path, folder)):
        os.makedirs(os.path.join(path, folder))
    path = os.path.join(path, folder)
    # fasta = open(fasta, newline='')
    fasta.seek(0)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta,'fasta'))
    records=[]
    for k,record in fasta_dict.items():
        newname=MOTUDict[k]
        record.id=newname+'_'+k
        record.description=''
        records.append(record)
    SeqIO.write(records,os.path.join(path, 'MOTU.fasta'),'fasta')

    newfasta2 = sorting_data(path, os.path.join(path,'MOTU.fasta'), outname='MOTU_sorted.fasta')  # Checking fasta
    if newfasta2.check_fasta() is False:
        sys.exit("Error: different length sequences, please check your alignment")
    newfasta2.sort_fasta()  # sorting fasta
    
def print_MOTU_X(path,MOTUs_df,k):
    """Print the delimitation results and write it to csv file"""
    MOTUprint = MOTUs_df.groupby(k).groups
    for k,v in MOTUprint.items():
        logging.info('#####'+k+'#####')
        logging.info('%s' % ', '.join(map(str, v)))   
        logging.info('')
    MOTUs_df.to_csv(os.path.join(path, 'MOTU_list.csv')) 
        
def extract_MOTU_nominal(basepath, fasta, gen=1, sp=2):
    if not os.path.exists(os.path.join(basepath,'Nominal/')):
        os.makedirs(os.path.join(basepath, 'Nominal/'))  
    with open(fasta, newline='') as handle:
        fasta_seqIO = SeqIO.parse(handle, 'fasta')
        names=[records.id for records in fasta_seqIO]
        spnames=['_'.join(i.split('_')[:2]) for i in names]
        datas = {"Nominal": spnames}
        MOTUs = pd.DataFrame(datas, index=names)                   
        return MOTUs

def dict_to_matrian(basepath,AllMOTUs,fasta,gen,sp,dis,cmd=False,diagnostic=False,n_ind=3):
    MOTUs_dict=AllMOTUs.to_dict()
    disdic={}
    for k,v in MOTUs_dict.items():
        out_folder = os.path.join(basepath, k)
        logging.info('\n#####################\n' + k + ' MOTUs\n#####################\n')
        # sys.argv = [sys.argv[0]]
        MOTU_X(basepath, fasta, out_folder, v, k)
        print_MOTU_X(basepath,AllMOTUs[[k]],k)
        logging.info('')
        distances = Matrian.main(os.path.join(basepath, out_folder), os.path.join(basepath,out_folder,'MOTU_sorted.fasta'), gen, sp, dis)
        disdic[k]=distances
        # if diagnostic == True:
        #     Diagnoser.main(os.path.join(basepath, out_folder), 'MOTU_sorted.fasta', n_ind)
        if cmd==True:
            distances.print_sum()
            distances.plot_max_min()
            distances.plot_freq() 
            if diagnostic==True:
                logging.info(("\n" + "### Summary of diagnostic character result ###\n"))
                diagnostics_cmd(basepath,k,n_ind)
    if len(MOTUs_dict)>1:
        return disdic
    else:
        return distances
    
def diagnostics(basepath,k,n_ind=3):
    out_folder= os.path.join(basepath, k)
    datas=Diagnoser.main(os.path.join(basepath, out_folder), 'MOTU_sorted.fasta', n_ind)
    return datas

def diagnostics_cmd(basepath,k,n_ind=3):
    datas=diagnostics(basepath,k,n_ind)
    logging.info(datas.summary)

class Compare:
    """Compare the delimitation resuls and produce Consensus MOTUs
    Parameters
    ----------
    path: str
        The path to folder with the input file.
    fasta: str
        The name of sorted fasta file.
    CompList: str
        list of delimitation methods used.       
    """
    def __init__(self, path, fasta, CompList):
        if not os.path.exists(os.path.join(path, 'Compare/')):
            os.makedirs(os.path.join(path, 'Compare/'))
        self.path = path
        self.CompList = CompList
        self.CompList = self.complist_read()
        self.compared = self.compare()
        self.fasta=fasta
        # self.fasta = open(fasta, newline='')
        
    def complist_read(self):
        if self.CompList == 'All':
            comp=[]
            directories=os.listdir(self.path)
            for i in directories:
                if os.path.exists(os.path.join(self.path, i, 'MOTU_list.csv')):
                    if i != 'Compare':
                        comp.append(i)
            if 'Nominal' in comp:
                comp.remove('Nominal')
                comp=['Nominal']+comp
            CompList = comp                    
        else:
            CompList = self.CompList.split(',')
        return CompList
            

    def compare(self):
        """Parse the different results and generate a file with Consensus Motus""" 
        with open(os.path.join(self.path, 'Compare','ALL_MOTUs.txt'), 'w+') as compare_file:
            AllMOTUs_df = pd.DataFrame() #dataframes with all motus
            AllMOTUs_dict={} #ditionary with all motus
            lists = [] #individuals by motu
            motus_name = [] #names of motus           
            n_analysis = 0
            for i in self.CompList:
                i_folder = i
                if os.path.exists(os.path.join(self.path, i_folder, 'MOTU_list.csv')):
                    n_analysis += 1
                    MOTUs_df = pd.read_csv(os.path.join(self.path, i_folder, 'MOTU_list.csv'), encoding="ISO-8859-1",index_col=0)
                    AllMOTUs_df[i] = MOTUs_df
                    MOTUs = MOTUs_df.groupby(i).groups
                    # print(MOTUs_df)
                    for k,v in MOTUs.items():
                        motus_name.append(k+'_('+i+')')   
                        lists.append(list(v))
                        AllMOTUs_dict[k+'_('+i+')']=list(v)
            used=[]
            sameind={}
            for k,v in AllMOTUs_dict.items():
                if k not in used:
                    concordants=[]
                    for l,w in AllMOTUs_dict.items():
                        if len(set(v)) == len(set(w)) and len(set(v)) == len(set.intersection(set(v), set(w))):
                            concordants.append(l)
                            used.append(l)
                    key='&'.join(concordants)
                    sameind[key]=v
            for k,v in sameind.items():                    
                compare_file.write(('%s' % ' & '.join(map(str, k.split('&')))) + '\n')
                compare_file.write(('%s' % ', '.join(map(str, v))) + '\n')

        num_motu = 1
        TC = {}
        MC = {}
        TD = {}
        MD = {}
        namedf=[]
        motudf=[]
        # print(sameind)
        if 'Nominal' in self.CompList and os.path.exists(os.path.join(self.path, 'Nominal','MOTU_list.csv')):
            n_analysis += -1
            for k,v in sameind.items():
                if '(Nominal)' in k:
                    # print('motus vistas con nominal',k)
                    n_used = len(k.strip().split("&"))-1
                    if n_used == n_analysis:
                        TC[k]=v
                    if n_used < n_analysis and n_used >= (n_analysis/float(2)):
                        MC[k]=v
                else:
                    n_used = len(k.strip().split("&"))
                    if n_used == n_analysis:
                        TD[k]=v
                    if n_used < n_analysis and n_used >= n_analysis / float(2):
                        MD[k]=v 
           
            logging.info('### All MOTUs match with the taxonomy ###\n')
            for k,v in TC.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n')
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)
                num_motu += 1 
                
            logging.info('### Most MOTUs agree with the taxonomy ###\n') 
            for k,v in MC.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n')
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)
                num_motu += 1 
                
            logging.info('### All MOTUs do not match the taxonomy###\n')
            for k,v in TD.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n')
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)
                num_motu += 1 
                
            logging.info('### Most MOTUs do not match the taxonomy ###\n')                    
            for k,v in MD.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n')
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)
                num_motu += 1 
                
            datas = {"Consensus": namedf}
            prints = pd.DataFrame(datas, index=motudf)
            prints=prints.sort_index(ascending=True)                     
            prints.to_csv( os.path.join(self.path, 'Compare','MOTU_list.csv'))

            
        else:
            for k,v in sameind.items():
                n_used = len(k.strip().split("&"))
                if n_used == n_analysis:
                    TC[k]=v
                if n_used < n_analysis and n_used >= (n_analysis/float(2)):
                    MC[k]=v
            logging.info('### MOTUs totally agree ###\n')
            for k,v in TC.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n') 
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)
                num_motu += 1 
                
            logging.info('### Most MOTUs match the same delimitation ###\n')            
            for k,v in MC.items():
                logging.info('Consensus MOTU ' + str(num_motu) + ' [' + k + ']')
                logging.info(', '.join(v)+'\n')
                namedf.extend(['MOTU_'+str(num_motu)]*len(v))
                motudf.extend(v)      
                num_motu += 1 
                
            datas = {"Consensus": namedf}
            prints = pd.DataFrame(datas, index=motudf)  
            prints=prints.sort_index(ascending=True)                  
            prints.to_csv( os.path.join(self.path, 'Compare','MOTU_list.csv'))
            
        AllMOTUs_df['Consensus']=prints #AGREGAR CONSENSUS LUEGO DEL WARNING
        
        # print(AllMOTUs_df)
        return AllMOTUs_df

    def Compare_warning(self):
        """Rise warnings in case of overlapping Consensus MOTUs and try to resolve it before continue the analysis""" 
        lists_w_temp = []
        lists_W = []
        motus_name_W = []
        motus_memb_W = []
        a = '########################################\n'
        with open(os.path.join(self.path, 'Compare','Summary_MOTUs.txt'), 'r+') as summary_file:
            summary_file.seek(0)
            for line in summary_file:
                if line.startswith('Consensus'):
                    motus_name_W.append(line.strip())
                    member = next(summary_file)
                    lists_w_temp.append(member.strip().split(', '))
                    motus_memb_W.append(member)
            for i in range(len(lists_w_temp)):  # MC vs MC
                for j in range(len(lists_w_temp)):
                    if i != j and i > j:
                        if set(lists_w_temp[i]) & set(lists_w_temp[j]) != set([]):
                            lists_W.append([i, j])

            if lists_W != []:
                logging.info('\n' + '### WARNING: overlapping MOTUS  ###\n')
                summary_file.write('\n' + a + 'WARNING: overlapping MOTUS\n' + a + '\n')
                n = 0
                for i in lists_W:
                    n += 1
                    logging.info('\n' + 'Overlapping consensus MOTUs ' + str(n) + '\n')
                    summary_file.write('\n' + 'Overlapping consensus MOTUs ' + str(n) + '\n')
                    logging.info(motus_name_W[i[0]].strip())
                    summary_file.write(motus_name_W[i[0]])
                    logging.info(motus_memb_W[i[0]])
                    summary_file.write(motus_memb_W[i[0]] + '\n')
                    logging.info(motus_name_W[i[1]].strip())
                    summary_file.write(motus_name_W[i[1]])
                    logging.info(motus_memb_W[i[1]])
                    summary_file.write(motus_memb_W[i[1]] + '\n')

                warn_chos = str(input(
                    'Do you want to choice between the overlapping MOTUs before calculate DNA barcoding statistics? (Y/N)').strip())
                while warn_chos != 'Y' and warn_chos != 'y' and warn_chos != 'n' and warn_chos != 'N':
                    logging.info('Error: Answer must be Y or N')
                    warn_chos = str(input(
                        'Do you want to chose between the overlapping MOTUs before calculate DNA barcoding statistics? (Y/N)').strip())
                if warn_chos == 'Y' or warn_chos == 'y':
                    rmvd, acpt, num_del = [], [], []
                    for i in lists_W:
                        motuA = motus_name_W[i[0]].strip()
                        motuB = motus_name_W[i[1]].strip()
                        text_cho = '1)' + motuA + ' or 2)' + motuB + ' (1/2)'
                        motu_cho = int(input(text_cho))
                        while motu_cho != 1 and motu_cho != 2:
                            logging.info('Error: Answer must be 1 or 2')
                            motu_cho = int(input(text_cho))
                        logging.info('\nMOTU chosen ' + motus_name_W[i[motu_cho - 1]] + '\n')
                        summary_file.write('MOTU chosen ' + motus_name_W[i[motu_cho - 1]] + '\n')
                        if motu_cho - 1 == 0:
                            motu_del = 1
                        else:
                            motu_del = 0
                        rmvd.append(motus_name_W[i[motu_del]])
                        acpt.append(motus_name_W[i[motu_cho - 1]])
                        num_del.append(motu_del)

                    while set(rmvd) & set(acpt) != set([]):
                        logging.info('You cant remove and acepted the same MOTU!!!!\n')

                        rmvd, acpt = [], []
                        for i in lists_W:
                            motuA = motus_name_W[i[0]].strip()
                            motuB = motus_name_W[i[1]].strip()
                            text_cho = '1)' + motuA + ' or 2)' + motuB + ' (1/2)'
                            motu_cho = int(input(text_cho))
                            while motu_cho != 1 and motu_cho != 2:
                                logging.info('Error: Answer must be 1 or 2')
                                motu_cho = int(input(text_cho))
                            logging.info('\nMOTU chosen ' + motus_name_W[i[motu_cho - 1]] + '\n')
                            summary_file.write('MOTU chosen ' + motus_name_W[i[motu_cho - 1]] + '\n')
                            if motu_cho - 1 == 0:
                                motu_del = 1
                            else:
                                motu_del = 0
                            rmvd.append(motus_name_W[i[motu_del]])
                            acpt.append(motus_name_W[i[motu_cho - 1]])

                    for ind, i in enumerate(lists_W):
                        with open(os.path.join(self.path, 'Compare','Compare_MOTU_list.txt'), 'r+') as out_list:
                            list_out = out_list.readlines()
                        with open(os.path.join(self.path, 'Compare','Compare_MOTU_list.txt'), 'w+') as out_list:
                            for ind2, line in enumerate(list_out):
                                if line.startswith('MOTU'):
                                    new_line = list_out[ind2 + 1]
                                    if new_line == motus_memb_W[i[num_del[ind]]]:
                                        pass
                                    else:
                                        out_list.write(line)
                                        out_list.write(new_line)

                elif warn_chos == 'N' or warn_chos == 'n':
                    logging.info('You need to choice a consensus MOTUS to calculate distances and graphs')
                    exit()

    def MOTU_renameFasta_Compare(self):
        """Rename the fasta name sequences using Consensus MOTUs""" 
        self.fasta.seek(0)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(self.fasta,'fasta'))
        records=[]
        for k,record in fasta_dict.items():
            MOTUs_dict=self.compared.to_dict()
            newname=MOTUs_dict['Consensus'][k]
            record.id=newname+'_'+k
            record.description=''
            records.append(record)
        SeqIO.write(records,os.path.join(self.path, 'Compare','MOTU.fasta'),'fasta')            


def plot_compare_tree(path, tree, motudf, nocons=False, names=True, save=False):
    """To plot the tree graph with MOTUs results
    Parameters
    ----------
    path: str
        The path to folder with the input file.       
    tree: str
        The tree file. 
    motudf: df
        DataFrame with all motus data
    names: boolean
        For hide names leaf in the tree, default=True to show the names.          
    """    
    tree = toytree.tree(tree)  # ,tree_format=1)
    if nocons == True:
        del motudf['Consensus']
    AllMOTUs_dict = motudf.to_dict()

    tree = tree.mod.node_scale_root_height(2 * len(AllMOTUs_dict))
    tips = (tree.get_tip_labels())
    canvas, axes, mark = tree.draw(
        width=800,
        height=30 * len(tips),
        tip_labels=False,  # hide labels
        tip_labels_align=True,
        tip_labels_style={"-toyplot-anchor-shift": "80px"},
        scalebar=True,
    );
    xsep = 0.5
    n = 0
    # add rectangles for delimitation
    for k in AllMOTUs_dict:
        n += 1
        coord = tree.get_tip_coordinates()
        coor_start = [coord[0][1]]
        coor_end = []
        for i, j in enumerate(tips):
            if i < len(tips) - 1:
                if AllMOTUs_dict[k].get(j) != AllMOTUs_dict[k].get(tips[i + 1]):
                    coor_end.append(coord[i][1])
                    coor_start.append(coord[i + 1][1])
        coor_end.append(coord[len(tips) - 1][1])
        for i in range(len(coor_start)):
            if k == 'Nominal':
                axes.rectangle(xsep * n, xsep * n + 0.2,
                               coor_start[i] - 0.35,
                               coor_end[i] + 0.35,
                               color='#0000ff',
                               );
            elif k == 'Consensus':
                axes.rectangle(xsep * n, xsep * n + 0.2,
                               coor_start[i] - 0.35,
                               coor_end[i] + 0.35,
                               color='#ff0000',
                               );
            else:
                axes.rectangle(xsep * n, xsep * n + 0.2,
                               coor_start[i] - 0.35,
                               coor_end[i] + 0.35,
                               color='#00000',
                               );
    # add tip labels
    ypos = range(len(tree))
    xpos = [xsep * n + 0.5] * len(tips)

    tipstyle = {"font-size": "12px",
                "text-anchor": "start",
                "fill": "#00000"}
    axes.text(xpos, ypos, tips,
              style=tipstyle,
              )
    # add analyses names
    labels = [k for k in AllMOTUs_dict]
    tipstyle = {"text-anchor": "start", "-toyplot-anchor-shift": "0"}
    for i in labels:
        axes.text(
            [(x + 1) * xsep + 0.1 for x in range(len(labels))],
            [len(tips)] * len(labels),
            labels,
            angle=90,
            color='#00000',
            style=tipstyle
        );
    if 'Nominal' in AllMOTUs_dict:
        num_tax = []
        num_contigous = []
        num_con = 0
        for v in tips:
            v=v.split('_')[:2]
            if v not in num_tax:
                num_tax.append(v)
            if v in num_contigous:
                pass
            else:
                num_contigous = [v]
                num_con += 1
        if len(num_tax) != num_con:
            logging.info('#####\nWarning: Nominal species not contigous in the tree.\n#####\n')
    if save==True:
        toyplot.svg.render(canvas, os.path.join(path, "Compare","compare_tree_lines.svg"))            
    # return(canvas)                   

    
def reading_data(fasta,tree=None,CODE=None):
    basepath=os.path.dirname(fasta)
    targets = logging.StreamHandler(sys.stdout), logging.FileHandler(
        os.path.join(basepath,'SPdel_output.log'))  # from https://stackoverflow.com/questions/24204898/python-output-on-both-console-and-file/24206109#24206109
    logging.basicConfig(format='%(message)s', level=logging.INFO, handlers=targets)
    logging.info('\n\n############################################################################\n')
    logging.info('SPdel v2.0 - Species delimitation and statistics for DNA Barcoding data sets\n')
    logging.info('############################################################################\n')    
    inputs = sorting_data(basepath, fasta, typeCODE=CODE,tree=tree)  # Checking fasta
    if inputs.check_fasta() is False:
        sys.exit("Error: different length sequences, please check your alignment")
    else:
        logging.info('Sequences are aligned (same size)')                
    inputs.describe_fasta()  # describing fasta
    inputs.sort_fasta()  # sorting fasta
    if tree != None:
        if inputs.check_tree() is False:
            sys.exit("Error: different names in tree, please check your tree and fasta alignment")
    return inputs

def run_nominal(basepath,inputs,gen=1,sp=2,dis='k'):
    if inputs.check_pattern() == False:
        sys.exit('Please rename your file using "genera_especies_individual" format or provide a nominal list file')
    nMOTUs=extract_MOTU_nominal(basepath, os.path.join(basepath,'Nominal/sorted.fasta'))
    distances=dict_to_matrian(basepath,nMOTUs,inputs.fasta,gen,sp,dis)
    return distances

def run_PTP(basepath,inputs,gen=1,sp=2,dis='k'):
    sys.argv = [sys.argv[0]]
    ptp(basepath, inputs.tree)
    PTP_df=MOTU_listPTP(os.path.join(basepath, 'PTP'),'PTP.PTPhSupportPartition.txt','PTP')
    distances=dict_to_matrian(basepath,PTP_df,inputs.fasta,gen,sp,dis)
    return distances
    
def run_bPTP(basepath,inputs,niter='10000',sample='100',burnin='0.1',gen=1,sp=2,dis='k'):
    sys.argv = [sys.argv[0]]
    bptp(basepath,inputs.tree,niter,sample,burnin)
    bPTP_df=MOTU_listPTP(os.path.join(basepath, 'bPTP'),'bPTP.PTPhSupportPartition.txt','bPTP')
    distances=dict_to_matrian(basepath,bPTP_df,inputs.fasta,gen,sp,dis)
    return distances
    
def run_GMYC(basepath,inputs,gen=1,sp=2,dis='k'):
    sys.argv = [sys.argv[0]]
    gmyc_run(basepath,inputs.tree)
    GMYC_df=MOTU_listGMYC(os.path.join(basepath, 'GMYC'),'GMYC_MOTU.txt','GMYC')
    distances=dict_to_matrian(basepath,GMYC_df,inputs.fasta,gen,sp,dis) 
    return distances

def run_PTPList(basepath,inputs,PTPList,gen=1,sp=2,dis='k'):
    PTP_df=MOTU_listPTP(basepath,PTPList,'PTP')
    distances=dict_to_matrian(basepath,PTP_df,inputs.fasta,gen,sp,dis)
    return distances
    
def run_bPTPList(basepath,inputs,bPTPList,gen=1,sp=2,dis='k'):
    bPTP_df=MOTU_listPTP(basepath,bPTPList,'bPTP')  
    distances=dict_to_matrian(basepath,bPTP_df,inputs.fasta,gen,sp,dis)
    return distances
    
def run_GMYCList(basepath,inputs,GMYCList,gen=1,sp=2,dis='k'):
    GMYC_df=MOTU_listGMYC(basepath,GMYCList,'GMYC')  
    distances=dict_to_matrian(basepath,GMYC_df,inputs.fasta,gen,sp,dis)
    return distances

def run_csvList(basepath,inputs,XList,gen=1,sp=2,dis='k',cmd=False,diagnostic=False,n_ind=3):
    if XList == None:
        sys.exit('Please provide a MOTUs List')
    AllMOTUs = pd.read_csv(XList, encoding="ISO-8859-1",index_col=0)
    check_file = check_MOTUlist(inputs.fasta, AllMOTUs)
    if check_file == False:
        sys.exit('Error: Different individuals in MOTUs List provided')
    distances=dict_to_matrian(basepath,AllMOTUs,inputs.fasta,gen,sp,dis,cmd,diagnostic,n_ind)
    return distances

def run_comparison(basepath,inputs,CompList,nocons=False,gen=1,sp=2,dis='k'):
    logging.info('\n#####################\n Consensus MOTUs\n#####################\n')
    comp_analize = Compare(basepath, inputs.fasta, CompList)
    # comp_analize.Compare_warning()
    comp_analize.MOTU_renameFasta_Compare()
    # if inputs.tree != None:
    #     plot_compare_tree(basepath, inputs.tree, comp_analize.compared,nocons)
    distances=Matrian.main(os.path.join(basepath, 'Compare/'), os.path.join(basepath, 'Compare','MOTU.fasta'), gen, sp, dis)
    return comp_analize.compared,distances
    # if 'd' in a:
    #     Diagnoser.main(os.path.join(basepath, 'Compare/'), 'Compare_MOTU.fasta', n_ind)          
    
def run(fasta,a,tree,CODE=None,dis='k',niter='10000',sample='100',burnin='0.1',gen=1,sp=2,n_ind=3,nocons=False,XList=None,PTPList=None,bPTPList=None,GMYCList=None,CompList=None,diagnostic=False):
    basepath=os.path.dirname(fasta)
    inputs=reading_data(fasta,tree,CODE)
    if CODE != None:
        inputs.stop_codon()   
    if 'n' in a:
        distances=run_nominal(basepath,inputs,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq()
        if diagnostic==True:
            datas=diagnostics_cmd(basepath,'Nominal',n_ind)
            logging.info(datas[0])
    if 'P' in a:
        distances=run_PTP(basepath,inputs,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq()
        if diagnostic==True:
            diagnostics_cmd(basepath,'PTP',n_ind)
    if 'T' in a:
        distances=run_bPTP(basepath,inputs,niter,sample,burnin,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq()
        if diagnostic==True:
            diagnostics_cmd(basepath,'bPTP',n_ind)
    if 'G' in a:
        distances=run_GMYC(basepath,inputs,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq()
        if diagnostic==True:
            diagnostics_cmd(basepath,'GMYC',n_ind)          
    if 'p' in a: #to do: check individual file
        distances=run_PTPList(basepath,inputs,PTPList,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq() 
        if diagnostic==True:
            diagnostics_cmd(basepath,'PTP',n_ind)
    if 't' in a: #to do: check individual file
        distances=run_bPTPList(basepath,inputs,bPTPList,gen,sp,dis)
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq()
        if diagnostic==True:
            diagnostics_cmd(basepath,'bPTP',n_ind)
    if 'g' in a: #to do: check individual file
        distances=run_GMYCList(basepath,inputs,GMYCList,gen,sp,dis)         
        distances.print_sum()
        distances.plot_max_min()
        distances.plot_freq() 
        if diagnostic==True:
            diagnostics_cmd(basepath,'GMYC',n_ind)
    if 'x' in a:
        distances=run_csvList(basepath,inputs,XList,gen,sp,dis,cmd=True,diagnostic=False,n_ind=3)
    if CompList != None:
        distances=run_comparison(basepath,inputs,CompList,nocons,gen,sp,dis)
        if inputs.tree != None:
            plot_compare_tree(basepath, inputs.tree, distances[0],nocons,save=True)            
        distances[1].print_sum()
        distances[1].plot_max_min()
        distances[1].plot_freq()     

def print_options():
    """Print avaliable options""" 
    print('')
    print('SPdel v1.0 - Species delimitation and statistics for DNA Barcoding data sets')
    print('')
    print(
        'The sequences name should be separate for "_" (e.g. Genus_species_individual) or use -N option for rename sequences')
    print('')
    print(
        'usage: ./SPdel.py path_to_files/ fasta_file -n -P PTP_File -t tree_file -X MOTUList1.txt,MOTUList2.txt -C p,MOTUList1,MOTUList2')
    print('usage: ./SPdel.py path_to_files/ fasta_file -n -P -G -T -t tree_file')
    print('usage: ./SPdel.py path_to_files/ fasta_file -n -distance p -code VER')
    print('usage: ./SPdel.py path_to_files/ fasta_file -P PTP_File -G GMYC_File -T bPTP_File -t tree_file -B BIN_file')
    print('usage: ./SPdel.py path_to_files/ fasta_file -t tree_file -C n,p,t,MOTUList1')
    print("Options:")
    print("    -n           For nominal analysis.\n")
    print("    -distance    Substitution model, k for K2p or p for p-distance (default=k)")
    print("    -t           Specify the name of the input newick tree for PTP, bPTP and GMYC analysis.\n")
    print("    -N           Specify the text file including the nominal names for rename the sequences.\n")
    print("    -P           Specify the PTP output file.\n")
    print("    -G           Specify the GMYC output file.\n")
    print("    -T           Specify the bPTP output file.\n")
    print("    -B           Specify the text file including the BIN names obtained from BOLD.\n")
    print("    -X           Specify the text file including the MOTUs names obtained from any external method.\n")
    print("    -nocons      Dont draw consensus in tree graphic.\n")
    print("    -D           For Diagnostic character analysis.\n")
    print(
        "    -C           Specify the type of analisys to be compared, include n for nominal, p for PTP, t for bPTP, b for BIN, and any filename used in X option for external MOTU list.\n")
    print("    -code        Specify the genetic code used to test stop codon, VER or INV.\n")
    print("Options for nominal: \n")
    print('    -gen         Position of the genus name in the sequence name when split by "_" (default=1).\n')
    print('    -sp          Position of the species name in the sequence name when split by "_" (default=2)\n')
    print("Options for bPTP: \n")
    print("    -n_iter       Number of iteration for bPTP analysis (default=10000)\n")
    print("    -sample      Number of sampling for bPTP analysis (default=100)\n")
    print("    -burnin      Burnin for bPTP analysis (default=0.1)\n")
    print("Options for diagnostic character: \n")
    print(
        "    -n_ind       Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3)\n")


def main(argu=None):
    """main function to run all analysis
    Parameters
    ----------
    path: str
        The path to folder with the input file.
    fasta: str
        The name of sorted fasta file.
    gen: int
        genus position on sequence names splited by "_".        
    sp: int
        genus position on sequence names splited by "_".         
    distance: str
        Substitution model, k for K2p or p for p-distance (default=k). 
    niter: int
        number of iterations for MCMC
    sample: int
        sampling for MCMC. 
    burnin: float
        burnin for MCMC.        
    BinList: str
        BIN species delimitation results file precalculated.   
    XList: str
        Any species delimitation results file precalculated.   
    NomList: str
        Nominal species list.           
    CompList: str
        species delimitation results file precalculated.   
    PTPList: str
        PTP species delimitation results file precalculated.  
    GMYCList: str
        GMYC species delimitation results file precalculated.
    bPTPList: str
        bPTP species delimitation results file precalculated.
    n_ind: int
        Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3).   
    CODE: str
        Type of genetic code used to check stop codon, use 'VER' or 'INV'           
    """
    #    orig_stdout = sys.stdout
    #    f=open(path+'output.txt','w')
    #    sys.stdout = f
    if len(sys.argv) < 3:
        print_options()
        sys.exit()
    if argu is not None:
        for i in argu:
            sys.argv.append(i)
    #    print sys.argv
    # path = sys.argv[1]
    fasta = sys.argv[1]
    a = ''
    gen = 1
    sp = 2
    tree = None
    dis = 'k'
    niter = '10000'
    sample = '100'
    burnin = '0.1'
    PTPList = None
    bPTPList=None
    GMYCList=None
    XList = None
    CompList = None
    n_ind = 3
    CODE = None
    nocons = False
    diagnostic=False

    for i in range(len(sys.argv)):
        if sys.argv[i] == '-n':
            a += 'n'
        elif sys.argv[i] == '-t':
            i = i + 1
            tree = sys.argv[i]
        elif sys.argv[i] == '-distance':
            i = i + 1
            dis = sys.argv[i]
        elif sys.argv[i] == '-D':
            diagnostic=True
        elif sys.argv[i] == '-P':
            j = i + 1
            if j == len(sys.argv):
                a += 'P'
            elif sys.argv[j].startswith('-'):
                a += 'P'
            else:
                a += 'p'
                PTPList = sys.argv[j]            
        elif sys.argv[i] == '-G':
            j = i + 1
            if j == len(sys.argv):
                a += 'G'
            elif sys.argv[j].startswith('-'):
                a += 'G'
            else:
                a += 'g'
                GMYCList = sys.argv[j] 
        elif sys.argv[i] == '-T':
            j = i + 1
            if j == len(sys.argv):
                a += 'T'
            elif sys.argv[j].startswith('-'):
                a += 'T'
            else:
                a += 't'
                bPTPList = sys.argv[j] 
        elif sys.argv[i] == '-X':
            i = i + 1
            a += 'x'
            XList = sys.argv[i]
        # elif sys.argv[i] == '-N':
        #     i = i + 1
        #     NomList = sys.argv[i]
        elif sys.argv[i] == '-n_iter':
            i = i + 1
            niter = sys.argv[i]
        elif sys.argv[i] == '-sample':
            i = i + 1
            sample = sys.argv[i]
        elif sys.argv[i] == '-burnin':
            i = i + 1
            burnin = sys.argv[i]
        elif sys.argv[i] == '-gen':
            i = i + 1
            gen = sys.argv[i]
        elif sys.argv[i] == '-sp':
            i = i + 1
            sp = sys.argv[i]
        elif sys.argv[i] == '-n_ind':
            i = i + 1
            n_ind = sys.argv[i]
        elif sys.argv[i] == '-code':
            i = i + 1
            a += 'e'
            CODE = sys.argv[i]
        elif sys.argv[i] == '-C':
            j = i + 1
            a += 'C'
            if j == len(sys.argv):
                CompList = 'All'
            else:
                CompList = sys.argv[j]
        elif sys.argv[i] == '-nocons':
            i = i + 1
            nocons=True

    run(fasta,a,tree,CODE,dis,niter,sample,burnin,gen,sp,n_ind,nocons,XList,PTPList,bPTPList,GMYCList,CompList,diagnostic)


    # sys.stdout = orig_stdout
    # f.close()


if __name__ == "__main__":
    #    basepath = os.path.dirname(os.path.abspath(__file__))
    main()
