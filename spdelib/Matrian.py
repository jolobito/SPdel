# -*- coding: utf-8 -*-
"""
Created on Sun Apr 01 21:44:07 2018

@author: Jorge L. Ramirez
"""
# import csv
import os
# from sys import platform as sys_pf

# if sys_pf == "darwin":
#     import matplotlib

#     matplotlib.use("Qt5Agg")

import pandas as pd
import numpy as np
from Bio import SeqIO
from math import log, sqrt
import logging
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pio.templates.default = "plotly_white"


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

    def __init__(self, path, fasta, gen, sp, distance,upper):
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
        self.summary = self.Summary_distances()
        self.tograph = self.summary[2:]

    # genera lista de nombres de especies

    def name_sp(self, a, b):
        L = []
        sp = self.fasta_names
        for ind in sp:
            ind = ind.split("_")
            L.append(ind[a - 1] + "_" + ind[b - 1])
        Lname = list(set(L))
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
        ind_number = len(self.fasta_names)
        for i in range(ind_number - 1):
            for j in range(i + 1, ind_number - 1):
                ind1.append(self.fasta_names[i])
                ind2.append(self.fasta_names[j])
                if distance == "k":
                    num = self.k2Pdistance(
                        self.fasta_dict[self.fasta_names[i]
                                        ].seq, self.fasta_dict[self.fasta_names[j]].seq
                    )
                else:
                    num = self.pdistance(
                        self.fasta_dict[self.fasta_names[i]
                                        ].seq, self.fasta_dict[self.fasta_names[j]].seq
                    )
                if num == -0.0:
                    num = 0.0
                value.append(round(num, 5))  # code by neguinha!

        datas = {"ind1": ind1, "ind2": ind2, "distance": value}
        summ = pd.DataFrame(datas)
        return summ

    # calculate minimum, mean and maximum intraspecific distances
    def min_media_max_intra(self):
        distras = []
        min_intra_dict = {}
        mean_intra_dict = {}
        max_intra_dict = {}
        for sp in self.Lname:
            data_tmp = self.data.loc[(self.data['ind1'].str.startswith(
                sp + '_')) & (self.data['ind2'].str.startswith(sp + '_'))]
            min_intra_dict[sp] = data_tmp.distance.min()
            mean_intra_dict[sp] = data_tmp.distance.mean()
            max_intra_dict[sp] = data_tmp.distance.max()
            distras.extend(list(data_tmp['distance']))
        return min_intra_dict, mean_intra_dict, max_intra_dict, distras

    # calculate minimum, mean and maximum interspecific distances and NN
    def min_media_max_inter(self):
        dister = []
        min_inter_dict = {}
        mean_inter_dict = {}
        max_inter_dict = {}
        NN_dict = {}
        for sp in self.Lname:
            data_tmp1 = self.data.loc[(self.data['ind1'].str.startswith(
                sp + '_')) & (~self.data['ind2'].str.startswith(sp + '_'))]
            data_tmp2 = self.data.loc[(~self.data['ind1'].str.startswith(
                sp + '_')) & (self.data['ind2'].str.startswith(sp + '_'))]
            data_tmps = [data_tmp1, data_tmp2]
            data_tmp = pd.concat(data_tmps)
            min_inter_dict[sp] = data_tmp.distance.min()
            mean_inter_dict[sp] = data_tmp.distance.mean()
            max_inter_dict[sp] = data_tmp.distance.max()
            kpos = []
            list1 = data_tmp[data_tmp.distance ==
                             data_tmp.distance.min()]['ind1'].unique()
            list2 = data_tmp[data_tmp.distance ==
                             data_tmp.distance.min()]['ind2'].unique()
            for i in list1:
                if not i.startswith(sp + '_'):
                    kpos.append(i)
            for i in list2:
                if not i.startswith(sp + '_'):
                    kpos.append(i)
            kpos = list(set(['_'.join(i.split('_')[0:2]) for i in kpos]))
            NN_dict[sp] = ' '.join(kpos)
        data_tmp3 = self.data
        for sp in self.Lname:
            data_tmp3 = data_tmp3.loc[(~data_tmp3['ind1'].str.startswith(sp)) | (
                ~data_tmp3['ind2'].str.startswith(sp))]
        dister.extend(list(data_tmp3['distance']))
        return min_inter_dict, mean_inter_dict, max_inter_dict, NN_dict, dister

    def plot_max_min(self):  # max vs min graph
        data = self.tograph[0]
        data_max = max([data["inter"].max(), data["intra2"].max()])

        fig = px.scatter(data, x="intra2", y="inter",
                         labels={
                             "intra2": "Maximum intraspecific",
                             "inter": "Minimum to NN",
                         },
                         hover_name="names",
                         )
        fig.update_traces(marker=dict(color="LightSkyBlue",
                                      size=10, line=dict(width=1, color="DarkSlateGrey")))
        fig.add_shape(type="line", x0=0, y0=0, x1=data_max + 1, y1=data_max + 1,
                      line=dict(color="Grey", width=1, dash="dash"),
                      layer='below')
        fig.update_layout(
            title_text="Maximum intraspecific vs Minimum to NN",
            title_x=0.5,
            autosize=False, width=700, height=700
        )
        fig.write_image(os.path.join(self.path, "min_max.pdf"))
        fig.show()

    def plot_freq(self):  # Barcoding gap graph
        tra, ter = self.tograph[1][-1], self.tograph[2][-1]
        newBins_tra = len(set(tra)) // 3
        if newBins_tra == 0:
            newBins_tra = 1
        bin_width = 0.5
        max_tra = max(tra)
        max_ter = max(ter)
        counts_tra, bins_tra = np.histogram(tra, bins=np.arange(
            min(tra), max(tra) + bin_width, bin_width))
        counts_ter, bins_ter = np.histogram(tra, bins=np.arange(
            min(tra), max(tra) + bin_width, bin_width))
        max_counts_tra = max(counts_tra)
        max_counts_ter = max(counts_ter)

        fig = go.Figure()
        fig.add_trace(go.Histogram(x=tra,
                                   name="Intragroup distance",
                                   bingroup=1
                                   ))
        fig.add_trace(go.Histogram(x=ter,
                                   name="Intergroup distance",
                                   bingroup=1
                                   ))

        fig.update_yaxes(title_text="# of taxon pairs")
        fig.update_layout(
            xaxis_range=(0, max(max_ter, max_tra)),
            yaxis_range=(0, max(max_counts_ter, max_counts_tra))
        )
        fig.update_layout(
            barmode='overlay',
            title_text="DNA Barcoding gap",
            title_x=0.5,
            height=350,
            hovermode="x")

        fig.update_traces(opacity=0.75, xbins_size=0.25)
        fig.write_image(os.path.join(self.path, "barcoding_gap.pdf"))
        fig.show()

    def plot_heatmap(self,upper=None):
        dfinv=self.data[['ind2', 'ind1', 'distance']].copy()
        dfinv.rename(columns = {'ind2':'ind1', 'ind1':'ind2'}, inplace = True)
        dftot=pd.concat([self.data,dfinv])
        heatDF=dftot.pivot(index='ind1', columns='ind2', values='distance')
        heatDF=heatDF.fillna(0)
        if len(heatDF) >= 100:
            sizebackground=len(heatDF)*10 
        elif len(heatDF) >= 50:
            sizebackground=len(heatDF)*20
        elif len(heatDF) < 50:
            sizebackground=len(heatDF)*30             
        listnames=list(heatDF.columns)
        nameused=[]
        sppos = []
        spinf= []
        for pos,i in enumerate(listnames):
          i=i.split('_')[0:2]
          if i not in nameused:
            sppos.append(pos)
            spinf.append(pos-1)
            nameused.append(i)
        spinf=spinf[1:]
        spinf.append(len(listnames)-1)
        
        if upper==None:
            fig = px.imshow(heatDF, color_continuous_scale="teal_r", width=sizebackground+20, height=sizebackground)
        else:            
            fig = px.imshow(heatDF, color_continuous_scale="teal_r", width=sizebackground+20, height=sizebackground, zmin=0, zmax=upper)
        fig.update(data=[{'hovertemplate':'Individual 1:%{x}<br>Individual 2: %{y}<br><b>Distance: %{z}</b><extra></extra>'}])
        fig.update_xaxes(title_text='', tickprefix = ' ',tickmode='linear')
        fig.update_yaxes(title_text='', ticksuffix = '  ',tickmode='linear')
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', xaxis={'side':'top','tickangle':-90}, 
                          font={'family':'Arial','size':10,'color':'rgb(42, 86, 116)'},
            title={
                'text': "Genetic pairwise distance matrix",
                'y':0.05,
                'x':0.5,'font':{'size':20}})        
        #y lines        
        for i in range(len(sppos)):
          fig.add_shape(type="line",
              x0=-1, y0=sppos[i]-0.3, x1=-1, y1=spinf[i]+0.3,
              line=dict(color='darksalmon',width=8)
          )
        
        #x lines
        for i in range(len(sppos)):
          fig.add_shape(type="line",
              x0=sppos[i]-0.3, y0=-1, x1=spinf[i]+0.3, y1=-1,
              line=dict(color='darksalmon',width=8)      
          )
        fig.write_image(os.path.join(self.path, "heatmap.pdf"))
        fig.show()
                
        
        
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
        # print(e)
        # print(self.Lname)
        datas2 = {"inter": e, "intra2": y, "names": self.Lname}
        df = pd.DataFrame(datas2, index=a)
        title = ["minimum", "mean", "maximum"]
        tra = []
        ter = []
        tra.append(min(list(intra[0].values())))
        tra.append(round(sum(list(intra[1].values())) / len(intra[1]), 5))
        tra.append(max(list(intra[2].values())))
        ter.append(min(list(inter[0].values())))
        ter.append(round(sum(list(inter[1].values())) / len(intra[1]), 5))
        ter.append(max(list(inter[2].values())))
        datas3 = {"intra": tra, "inter": ter}
        tab = pd.DataFrame(datas3, index=title)

        return summ, tab, df, intra, inter

    def print_sum(self):
        logging.info(
            "Summary table (Name, mean intra, max intra, NN, distance to NN) in percentage"
        )
        logging.info(self.summary[0])
        logging.info("")
        logging.info(self.summary[1].transpose())
        return self.summary[0], self.summary[1].transpose()

    def print_summary(self):
        return self.summary[0]

    def print_summary_all(self):
        return self.summary[1].transpose()



def main(path, fasta_file, gen, sp, distance, upper=None, out_name=None, n=False):
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
    tmp = Matrian(path, fasta_file, gen, sp, distance,upper)
   # tmp.analyze()
    return tmp
    # tmp.graphics(tograph[0],tograph[1],tograph[2])

    # if path is not None:
    #     tmp.save_csv(fasta_file)


if __name__ == "__main__":
    main()
# import time
# h=time.time()
# tmp=Matrian('C:/Users/JORGE/OneDrive/Python/SPdel/SPdel_2023_lite/data',"C:/Users/JORGE/OneDrive/Python/SPdel/SPdel_2023_lite/data/Laemolyta/LA_nominal.fasta",1,2,'k')
# tmp.analyze()
# g=time.time()
# print(g-h)
# tmp.save_csv('C:/Users/ramir/OneDrive/Python/SPdel/Example//Schizodon/SC_COI_ALINHADO.fasta')
