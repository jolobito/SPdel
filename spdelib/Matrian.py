# -*- coding: utf-8 -*-
"""
Created on Sun Apr 01 21:44:07 2018

@author: Jorge L. Ramirez
"""
# import csv
import os
from sys import platform as sys_pf

if sys_pf == "darwin":
    import matplotlib

    matplotlib.use("Qt5Agg")

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import SeqIO
from math import log, sqrt
import logging
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


class Matrian:
    """calculate and print main genetic distances results.
    Parameters
    ----------
    path: str
        The path to folder to save outputs files.
    fasta: str
        The name of fasta file.
    gen: int
        genus position on sequence names splited by "_".
    sp: int
        genus position on sequence names splited by "_".
    distance: str
        Substitution model, k for K2p or p for p-distance (default=k).
    """

    def __init__(self, path, fasta, gen, sp, distance):
        self.path = path
        self.fasta = open(fasta, newline="")     
        self.fasta_seqIO = SeqIO.parse(self.fasta, "fasta")  # sequence object           
        self.fasta_names = [
            seq.id for seq in self.fasta_seqIO
        ]  # list of samples names
        self.fasta_dict = SeqIO.to_dict(
            SeqIO.parse(fasta, "fasta")
        )  # dictionary of sequences
        self.Lname = self.name_sp(gen, sp)  # list of species
        self.data = self.matrix(distance)  # matrix genetic distances
        self.tmp = self.Summary_distances()
        self.tograph = self.tmp[2:]        


    # genera lista de nombres de especies
    def name_sp(self, a, b):
        L = []
        sp = self.fasta_names
        for ind in sp:
            ind = ind.split("_")
            L.append(ind[a - 1] + "_" + ind[b - 1])
        Lname=list(set(L))
        Lname.sort()
        return Lname

    # adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
    def k2Pdistance(self, seq1, seq2):
        pairs = []
        for x in zip(seq1, seq2):
            if "-" not in x and "N" not in x:
                pairs.append(x)
        ts_count = 0
        tv_count = 0
        length = len(pairs)
        transitions = ["AG", "GA", "CT", "TC"]
        transversions = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]
        for (x, y) in pairs:
            if x + y in transitions:
                ts_count += 1
            elif x + y in transversions:
                tv_count += 1
        p = float(ts_count) / length
        q = float(tv_count) / length
        try:
            d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q)) * 100
        except ValueError:
            logging.info("Tried to take log of a negative number")
            return None
        return d

    # adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
    @staticmethod
    def pdistance(seq1, seq2):
        p = 0
        pairs = []
        for x in zip(seq1, seq2):
            if "-" not in x and "N" not in x:
                pairs.append(x)
        for (x, y) in pairs:
            if x != y:
                p += 1
        length = len(pairs)
        return (float(p) / length) * 100

    def matrix(self, distance):  # generate all genetic distance
        if distance == "k":
            logging.info("Using k2p distance\n")
        else:
            logging.info("Using p-distance\n")                   
        ind1, ind2, value, = [], [], []                    
        ind_number=len(self.fasta_names)
        for i in range(ind_number-1):
            for j in range(i+1,ind_number-1):
                ind1.append(self.fasta_names[i])
                ind2.append(self.fasta_names[j])
                if distance == "k":
                    num = self.k2Pdistance(
                        self.fasta_dict[self.fasta_names[i]].seq, self.fasta_dict[self.fasta_names[j]].seq
                    )                    
                else:               
                    num = self.pdistance(
                        self.fasta_dict[self.fasta_names[i]].seq, self.fasta_dict[self.fasta_names[j]].seq
                    )
                if num == -0.0:
                    num = 0.0               
                value.append(round(num, 5))  # code by neguinha!                
                
        datas = {"ind1": ind1, "ind2": ind2, "distance": value}
        summ = pd.DataFrame(datas)
        return summ
           
    def min_media_max_intra(self):  # calculate minimum, mean and maximum intraspecific distances
        distras=[]
        min_intra_dict = {}
        mean_intra_dict = {}
        max_intra_dict = {}
        for sp in self.Lname:
            data_tmp=self.data.loc[(self.data['ind1'].str.startswith(sp+'_')) & (self.data['ind2'].str.startswith(sp+'_'))]
            min_intra_dict[sp] = data_tmp.distance.min()
            mean_intra_dict[sp] = data_tmp.distance.mean()
            max_intra_dict[sp] = data_tmp.distance.max()
            distras.extend(list(data_tmp['distance']))
        return min_intra_dict, mean_intra_dict, max_intra_dict, distras

    def min_media_max_inter(self):  # calculate minimum, mean and maximum interspecific distances and NN
        dister=[]
        min_inter_dict = {}
        mean_inter_dict = {}       
        max_inter_dict = {}        
        NN_dict = {}
        for sp in self.Lname:
            data_tmp1=self.data.loc[(self.data['ind1'].str.startswith(sp+'_')) & (~self.data['ind2'].str.startswith(sp+'_'))]
            data_tmp2=self.data.loc[(~self.data['ind1'].str.startswith(sp+'_')) & (self.data['ind2'].str.startswith(sp+'_'))]
            data_tmps=[data_tmp1,data_tmp2]
            data_tmp=pd.concat(data_tmps)
            min_inter_dict[sp] = data_tmp.distance.min()        
            mean_inter_dict[sp] = data_tmp.distance.mean()
            max_inter_dict[sp] = data_tmp.distance.max()            
            kpos=[]
            list1=data_tmp[data_tmp.distance == data_tmp.distance.min()]['ind1'].unique()
            list2=data_tmp[data_tmp.distance == data_tmp.distance.min()]['ind2'].unique()
            for i in list1:
                if not i.startswith(sp+'_'):
                    kpos.append(i)
            for i in list2:
                if not i.startswith(sp+'_'):
                    kpos.append(i)                
            kpos = list(set(['_'.join(i.split('_')[0:2]) for i in kpos]))    
            NN_dict[sp] = ' '.join(kpos)           
        data_tmp3=self.data
        for sp in self.Lname:
            data_tmp3=data_tmp3.loc[(~data_tmp3['ind1'].str.startswith(sp)) | (~data_tmp3['ind2'].str.startswith(sp))]
        dister.extend(list(data_tmp3['distance']))
        return min_inter_dict, mean_inter_dict, max_inter_dict, NN_dict, dister
    


    def plot_max_min(self):  # max vs min graph
        sns.lmplot(x="intra2", y="inter", data=self.tograph[0], fit_reg=False)
        plt.title("Maximum intraspecific vs Minimum to NN")
        plt.xlabel("Maximum intraspecific")
        plt.ylabel("Minimum to NN")
        z = [self.tograph[0]["inter"].max(), self.tograph[0]["intra2"].max()]
        plt.axis([0, max(z) + 1, 0, max(z) + 1])
        lims = [0, max(z) + 1]
        plt.plot(lims, lims, ":k")
        plt.savefig(os.path.join(self.path, "min_max.pdf"))
        plt.show()
        # plt.clf()

    def plot_freq(self):  # Barcoding gap graph
        tra,ter=self.tograph[1][-1],self.tograph[2][-1]
        newBins_tra = len(set(tra)) // 3
        if newBins_tra == 0:
            newBins_tra = 1
        newBins_ter = len(set(ter)) // 3
        f, ax = plt.subplots(3, 1, sharex="col", sharey="all")
        sns.histplot(
            ter,
            bins=newBins_ter,
            color="b",
            kde=False,
            label="Intergruop distance",
            ax=ax[1],
        )
        sns.histplot(
            tra,
            bins=newBins_tra,
            color="r",
            kde=False,
            label="Intragroup distance",
            ax=ax[0],
        )
        sns.histplot(
            tra,
            bins=newBins_tra,
            color="r",
            kde=False,
            label="Intragroup distance",
            ax=ax[2],
        )
        sns.histplot(
            ter,
            bins=newBins_ter,
            color="b",
            kde=False,
            label="Intergruop distance",
            ax=ax[2],
        )
        ax[0].set_title("DNA Barcoding gap")
        ax[0].set_ylabel("# of taxon pairs")
        ax[1].set_ylabel("# of taxon pairs")
        ax[2].set_ylabel("# of taxon pairs")
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        ax[2].set_xlabel("Genetic Distance")
        f.savefig(os.path.join(self.path, "barcoding_gap.pdf"))
        plt.show()
        # plt.close(f)
        # plt.clf()     
        
    def Summary_distances(self):  # hace analisis en bloque
        a, b, c, d, e, y = [], [], [], [], [], []
        inter = self.min_media_max_inter()
        intra = self.min_media_max_intra()
        for sp in self.Lname:
            a.append(sp)
            b.append(intra[1].get(sp))
            c.append(intra[2].get(sp))
            d.append(inter[3].get(sp))
            e.append(inter[0].get(sp))
            if intra[2].get(sp) != intra[2].get(sp):
                y.append(0.0)
            else:
                y.append(intra[2].get(sp))
        datas = {"Mean": b, "Max": c, "NN": d, "DtoNN": e}
        summ = pd.DataFrame(datas, index=a)        
        ####data for plot max vc min graph####
        datas2 = {"inter": e, "intra2": y}
        df = pd.DataFrame(datas2, index=a)
        title = ["minimum", "mean", "maximum"]
        tra = []
        ter = []        
        tra.append(min(list(intra[0].values())))
        tra.append(round(sum(list(intra[1].values()))/len(intra[1]),5))
        tra.append(max(list(intra[2].values())))
        ter.append(min(list(inter[0].values())))
        ter.append(round(sum(list(inter[1].values()))/len(intra[1]),5))
        ter.append(max(list(inter[2].values())))        
        datas3 = {"intra": tra, "inter": ter}
        tab = pd.DataFrame(datas3, index=title)       

        return summ,tab,df,intra,inter
        
    def print_sum(self):
        logging.info(
            "Summary table (Name, mean intra, max intra, NN, distance to NN) in percentage"
        )        
        logging.info(self.tmp[0])
        logging.info("")
        logging.info(self.tmp[1].transpose())
        return self.tmp[0],self.tmp[1].transpose()
    
    def print_summary(self):
        return self.tmp[0]

    def print_summary_all(self):
        return self.tmp[1].transpose()   
    
    # def save_csv(self, out_name):
    #     out_name_csv = os.path.splitext(out_name)[0] + ".csv"
    #     with open(out_name_csv, "w", newline="") as out_file:
    #         writer = csv.writer(out_file, delimiter=",")
    #         writer.writerow(["name"] + self.fasta_names)
    #         for r in self.fasta_names:
    #             writer.writerow(
    #                 [r] + [value for value in self.data[r].values()]
    #             )


def main(path, fasta_file, gen, sp, distance, out_name=None, n=False):
    """A function to calculate and print main genetic distances results.
    Parameters
    ----------
    path: str
        The path to output folder.
    fasta_file: str
        The name of fasta file.
    gen: int
        genus position on sequence names splited by "_".
    sp: int
        genus position on sequence names splited by "_".
    distance: str
        Substitution model, k for K2p or p for p-distance (default=k).
    out_name: str
        The output for norminal list file.
    n: boolean
        True if nominal analysis.
    """
    tmp = Matrian(path, fasta_file, gen, sp, distance)
   #tmp.analyze()   
    return tmp
    # tmp.graphics(tograph[0],tograph[1],tograph[2])

    # if path is not None:
    #     tmp.save_csv(fasta_file)


if __name__ == "__main__":
    main()
# import time
# h=time.time()
# tmp=Matrian('C:/Users/jramirez/OneDrive/Python/SPdel/SPdel_github/data',"C:/Users/jramirez/OneDrive/Python/SPdel/SPdel_github/data/LA_nominal.fasta",1,2,'k')
# tmp.analyze()
# g=time.time()
# print(g-h)
# tmp.save_csv('C:/Users/ramir/OneDrive/Python/SPdel/Example//Schizodon/SC_COI_ALINHADO.fasta')
