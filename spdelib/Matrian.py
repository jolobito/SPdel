# -*- coding: utf-8 -*-
"""
Created on Sun Apr 01 21:44:07 2018

@author: jolob
"""
import csv
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
from collections import defaultdict
import logging


class Matrian:
    """calculate and print main genetic distances results.
    Parameters
    ----------
    path: str
        The path to folder with the input file.
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

    def matrix(self, distance):  # genera la matriz de distancia
        if distance == "k":
            logging.info("Using k2p distance\n")
        else:
            logging.info("Using p-distance\n")
        dmatrix = defaultdict(dict)
        for i in self.fasta_names:
            for j in self.fasta_names:
                if distance == "k":
                    num = self.k2Pdistance(
                        self.fasta_dict[i].seq, self.fasta_dict[j].seq
                    )
                    if num == -0.0:
                        num = 0.0
                    dmatrix[i][j] = round(num, 5)  # code by neguinha!
                if distance == "p":
                    num = self.pdistance(
                        self.fasta_dict[i].seq, self.fasta_dict[j].seq
                    )
                    if num == -0.0:
                        num = 0.0
                    dmatrix[i][j] = round(num, 5)  # code by neguinha!
        return dmatrix

    # genera lista de nombres de especies
    def name_sp(self, a, b):
        L = []
        sp = self.fasta_names
        for ind in sp:
            ind = ind.split("_")
            L.append(ind[a - 1] + "_" + ind[b - 1])
        Lname = []
        [Lname.append(key) for key in L if key not in Lname]
        Lname.sort()
        return Lname

    def mean_intra(self):  # da media intra de todos los grupos
        mean_intra_dict = {}
        for sp in self.Lname:
            ind_sp = [ind for ind in self.fasta_names if ind.startswith(sp)]
            L = []
            for ind in ind_sp:
                ind_dic = self.data[ind]
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if key.startswith(sp) and key is not ind
                ]
            if L == []:
                v = None
            else:
                v = sum(L) / len(L)
            mean_intra_dict[sp] = v
        return mean_intra_dict

    def max_intra(self):  # da maxima dentro de un grupo
        max_intra_dict = {}
        for sp in self.Lname:
            ind_sp = [ind for ind in self.fasta_names if ind.startswith(sp)]
            L = []
            for ind in ind_sp:
                ind_dic = self.data[ind]
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if key.startswith(sp)
                ]
            v = max(L)
            max_intra_dict[sp] = v
        return max_intra_dict

    def min_inter(self):  # da minima inter de um grupo
        min_inter_dict = {}
        NN_dict = {}
        for sp in self.Lname:
            ind_sp = [
                ind for ind in self.fasta_names if ind.startswith(sp)
            ]  # list of individuals of a given sp
            L = []
            K = []
            for ind in ind_sp:
                ind_dic = self.data[
                    ind
                ]  # dictionary of individuals genetic distances
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if not key.startswith(sp)
                ]
                K = K + [
                    key for key in ind_dic.keys() if not key.startswith(sp)
                ]
            v = min(L)
            kpos = L.index(v)
            min_inter_dict[sp] = v
            NN_dict[sp] = K[kpos]
        return min_inter_dict, NN_dict

    def min_media_max_intra(self):  # newused #review!! should be by group
        # mmm_intra_dict={}
        L = []
        for sp in self.Lname:
            ind_sp = [ind for ind in self.fasta_names if ind.startswith(sp)]
            for ind in ind_sp:
                ind_dic = self.data[ind]
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if key.startswith(sp) and key is not ind
                ]
        mintra = min(L)
        mediatra = round((sum(L) / len(L)), 5)
        maxtra = max(L)
        L.sort()
        return mintra, mediatra, maxtra, L

    def min_media_max_inter(self):  # newused #review!! should be by group
        L = []
        for sp in self.Lname:
            ind_sp = [
                ind for ind in self.fasta_names if ind.startswith(sp)
            ]  # list of individuals of a given sp
            for ind in ind_sp:
                ind_dic = self.data[
                    ind
                ]  # dictionary of individuals genetic distances
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if not key.startswith(sp)
                ]
        minter = min(L)
        mediater = round((sum(L) / len(L)), 5)
        maxter = max(L)
        L.sort()
        return minter, mediater, maxter, L

    def plot_max_min(self, df):  # max vs min graph
        sns.lmplot(x="intra2", y="inter", data=df, fit_reg=False)
        plt.title("Maximum intraspecific vs Minimum to NN")
        plt.xlabel("Maximum intraspecific")
        plt.ylabel("Minimum to NN")
        z = [df["inter"].max(), df["intra2"].max()]
        plt.axis([0, max(z) + 1, 0, max(z) + 1])
        lims = [0, max(z) + 1]
        plt.plot(lims, lims, ":k")
        plt.savefig(os.path.join(self.path, "min_max.pdf"))
        # plt.show()
        plt.clf()

    def plot_freq(self, tra, ter):  # Barcoding gap graph
        newBins_tra = len(set(tra)) // 3
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
        # plt.show(f)
        # plt.close(f)
        plt.clf()

    def analyze(self):  # hace analisis en bloque
        logging.info(
            "Summary table (Name, mean intra, max intra, NN, distance to NN) in percentage"
        )
        a, b, c, d, e, y = [], [], [], [], [], []
        mininter = self.min_inter()
        maxintra = self.max_intra()
        meanintra = self.mean_intra()
        for sp in self.Lname:
            a.append(sp)
            b.append(meanintra.get(sp))
            c.append(maxintra.get(sp))
            d.append(mininter[1].get(sp))
            e.append(mininter[0].get(sp))
            if maxintra.get(sp) == None:
                y.append(0.0)
            else:
                y.append(maxintra.get(sp))
        datas = {"Mean": b, "Max": c, "NN": d, "DtoNN": e}
        summ = pd.DataFrame(datas, index=a)
        # print(summ)
        logging.info(summ)
        datas2 = {"inter": e, "intra2": y}

        df = pd.DataFrame(datas2, index=a)

        logging.info("")
        tra = []
        ter = []
        title = ["minimum", "mean", "maximum"]
        L = self.min_media_max_intra()
        for i in range(len(L) - 1):
            tra.append(L[i])
        M = self.min_media_max_inter()
        for i in range(len(M) - 1):
            ter.append(M[i])
        datas3 = {"intra": tra, "inter": ter}
        tab = pd.DataFrame(datas3, index=title)
        # print(tab.transpose())
        logging.info(tab.transpose())

        ####Plot max vc min graph####

        self.plot_max_min(df)

        ####Plot frequencies graph####

        self.plot_freq(L[-1], M[-1])

    def save_csv(self, out_name):
        out_name_csv = os.path.splitext(out_name)[0] + ".csv"
        with open(out_name_csv, "w", newline="") as out_file:
            writer = csv.writer(out_file, delimiter=",")
            writer.writerow(["name"] + self.fasta_names)
            for r in self.fasta_names:
                writer.writerow(
                    [r] + [value for value in self.data[r].values()]
                )


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
    tmp.analyze()

    if path is not None:
        tmp.save_csv(fasta_file)

    if n == True:
        List = tmp.name_sp(gen, sp)
        with open(out_name, "w+") as nominal_file:
            for i in List:
                nominal_file.write(i + "\n")


if __name__ == "__main__":
    main()

# tmp=matrian('C:/Users/ramir/OneDrive/Python/SPdel/Example//Schizodon/',"SC_COI_ALINHADO.fasta",1,2,'k')
# tmp.analyze()
# tmp.save_csv('C:/Users/ramir/OneDrive/Python/SPdel/Example//Schizodon/SC_COI_ALINHADO.fasta')
