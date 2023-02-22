# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 00:29:53 2018

@author: jolob
"""

# import sys
# import plotly.plotly as py
# import sys
import os
import shutil
from math import ceil

import pandas as pd
import plotly

# import plotly.plotly as py
import plotly.graph_objs as go
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# import logging
# from Bio.Alphabet import IUPAC, Gapped
# alphabet = Gapped(IUPAC.ambiguous_dna)


class diagnostic_character:
    """all values necesary for diagnostic character plots
    
    Parameters
    ----------
    basepath: str
        The folder for outputs.    
    fasta: str
        The fasta file.
    """
    def __init__(self, basepath, fasta, n_ind=3, gen=1, sp=2):
        self.path = basepath
        self.fasta = fasta
        self.n_ind=int(n_ind)
        self.gen=gen
        self.sp=sp
        self.list_n_pos = self.snper()
        self.n_spec = self.specier(os.path.join(self.path, "tmp/"), "snps_all.fasta", self.gen - 1, self.sp - 1)
        self.spec_list = self.n_spec[1]
        self.n_ind_spec = self.n_spec[2]
        self.n_spec = self.n_spec[0]
        self.list_all = self.diagnoser()     
        self.summary_tmp = self.summary_all()
        self.summary = self.summary_tmp[0]
        self.spec_list2 = self.summary_tmp[1]

    # list_n_pos = snper(path, fasta)
    # n_spec = specier(os.path.join(path, "tmp/"), "snps_all.fasta", gen - 1, sp - 1)
    # spec_list = n_spec[1]
    # n_ind_spec = n_spec[2]
    # n_spec = n_spec[0]
    # list_all = diagnoser(path, n_spec, n_ind)        

    def snper(self):  # code modified from Jufisawa blog
        """calculate the positions with polymorphism.
        Parameters
        ----------
        path: str
            The path to folder with the input file.
        fasta: str
            The name of fasta file.
        Returns
        -------
        A list of polymorphic positions.
        """    
        if not os.path.exists(os.path.join(self.path, "tmp/")):
            os.makedirs(os.path.join(self.path, "tmp"))
        s_bases = ("A", "C", "G", "T")
        d_bases = ("R", "Y", "S", "W", "K", "M", "N")
        bases = s_bases + d_bases
        alig = AlignIO.read(os.path.join(self.path, self.fasta), "fasta")
        list_n_pos = []
        # snp = [SeqRecord(Seq('', s.seq.alphabet), id=s.id, description=s.description) for s in alig]
        snp = [
            SeqRecord(Seq(""), id=s.id, description=s.description) for s in alig
        ]
        snp = MultipleSeqAlignment(snp)
        for i in range(alig.get_alignment_length()):
            col = alig[:, i]
            col = col.upper()
            col = col.replace("-", "N")
            base_count = {b: col.count(b) for b in bases}
            genotypes = sorted(list(base_count.items()), key=lambda x: -x[1])
            genotypes = [b for b in genotypes if b[1] > 0]
            if len(genotypes) > 1:
                list_n_pos.append(i + 1)
                snp = snp + alig[:, i : i + 1]
        with open(os.path.join(self.path, "tmp/snps_all.fasta"), "w+") as snpfile:
            snpfile.write(format(snp,"fasta"))
        #    print list_n_pos
        return list_n_pos

    def specier(self,path,fasta,gen,sp):
        """return a tuple with basic information.
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
        Returns
        -------
        The number of species, a list with name of species, and a list with number of individual for each species.        
        """    
        alig = AlignIO.read(os.path.join(path, fasta), "fasta")
        spec = []
        n_spec = 0
        n_ind_spec = []
        for seq in alig:
            new_seq_id = seq.id.split("_")
            new_seq_id = new_seq_id[gen] + "_" + new_seq_id[sp]
            if new_seq_id not in spec:
                spec.append(new_seq_id)
                n_spec += 1
        for i in range(n_spec):
            n = 0
            with open(os.path.join(path, "spec" + str(i) + ".fasta"), "w+") as fasta:
                for seq in alig:
                    if seq.id.startswith(spec[i]):
                        fasta.write(">" + seq.description + "\n")
                        fasta.write(str(seq.seq).replace("-", "N") + "\n")
                        n += 1
            n_ind_spec.append(n)
    
        return n_spec, spec, n_ind_spec

    def diagnoser(self):
        """return a tuple with character diagnostic information.
        Parameters
        ----------
        path: str
            The path to folder with the input file.
        n_spec: str
            number of individuals for each species.    
        n_ind: int
            minimum number of individuals for species.
        Returns
        -------
        A list with a list for each polymorfic position with type of character of each species, and a list of list of snps for each position and each species.        
        """      
        lists_var = []
        for i in range(self.n_spec):
            lists_var.append(
                self.consensus_fasta(os.path.join(self.path, "tmp/"), "spec" + str(i) + ".fasta", self.n_ind)
            )
        len_lists_var=max([len(i) for i in lists_var if i is not None])
        list_all = []
        for j in range(len_lists_var):
            list_D = []
            N = 0
            n_var = []
            N_only = []
            for l in range(self.n_spec):
                if lists_var[l] is not None:
                    n_var.append(len(lists_var[l][j]))
                    for nucl in range(len(lists_var[l][j])):
                        if set(lists_var[l][j][nucl]) & set(["N"]) == set(["N"]):
                            N = 1
                    if lists_var[l][j] == ["N"]:
                        N_only.append(1)
                    else:
                        N_only.append(0)
            for i in range(self.n_spec):  # 10
                if lists_var[i] is not None:
                    list_tmp = []
                    #                    print lists_var[i][j]
                    for k in range(self.n_spec):  # 10
                        if lists_var[k] is not None:
                            if set(lists_var[i][j]) & set(lists_var[k][j]) != set(
                                []
                            ):
                                list_tmp.append(1)
                            else:
                                list_tmp.append(0)
                    count = 0
                    part = 0
                    n_only = 0
                    for b in list_tmp:
                        if b == 1:
                            count += 1
                    for c in N_only:
                        if c == 1:
                            n_only += 1
                    for d in range(len(list_tmp)):
                        if list_tmp[d] == 1 and n_var[d] == 1:
                            part += 1
                    if (
                        count == 1
                        and N == 0
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                    ):
                        list_D.append(1)  # 1 Diagnostic
                    elif (
                        count == 1
                        and N == 1
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                    ):
                        list_D.append(2)  # 2 Diagnostic or partial
                    elif (
                        count > 1
                        and N == 0
                        and part == 1
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                    ):
                        list_D.append(3)  # 3 Partial
                    elif (
                        count > 1
                        and N == 1
                        and part == 1
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                    ):
                        list_D.append(4)  # 4 Partial or uniformative
                    elif (
                        count == 1
                        and N == 1
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                        and n_only >= 1
                    ):
                        list_D.append(5)  # 5 Invalid
                    elif (
                        count > 1
                        and N == 1
                        and part == 1
                        and len(lists_var[i][j]) == 1
                        and lists_var[i][j] != ["N"]
                        and n_only >= 1
                    ):
                        list_D.append(5)  # 5 Invalid
                    else:
                        list_D.append(0)  # Uniniformative
            #            print list_D
            list_all.append(list_D)
        #    print lists_var
        return list_all, lists_var

    def consensus_fasta(self,path, fasta, n_ind):  # modified from biopython  #add bases degeneradas
        alig = AlignIO.read(os.path.join(path, fasta), "fasta")
        #    print alig
        consensus = []
        con_len = alig.get_alignment_length()
        n_seq = len(alig)
        if n_seq >= n_ind:
            for n in range(con_len):
                pos_bases = []
                for record in alig:
                    if n < len(record.seq):
                        if record.seq[n] not in pos_bases:
                            pos_bases.append(record.seq[n])
                consensus.append(pos_bases)
            #        print consensus
            return consensus
        else:
            #        print None
            return None

    def summary_all(self):
        """print a summary of kind of diagnostic character by each species. 
        Returns
        -------
        A list with names of species included in the analysis.        
        """        
        list_all=self.list_all[0]
        with open(os.path.join(self.path, "Diagnostic_character.txt"), "w+") as diagnostic:
            u, w, x, y, z = [], [], [], [], []
            spec_list2 = []
            n_ind_spec2 = []
            for i in range(len(self.spec_list)):
                if self.n_ind_spec[i] >= self.n_ind:
                    spec_list2.append(self.spec_list[i])
                    n_ind_spec2.append(self.n_ind_spec[i])
            for k in range(len(spec_list2)):
                u_sum, w_sum, x_sum, y_sum, z_sum = 0, 0, 0, 0, 0
                for j in list_all:
                    if j[k] == 1:
                        u_sum += 1
                    elif j[k] == 2:
                        w_sum += 1
                    elif j[k] == 3:
                        x_sum += 1
                    elif j[k] == 4:
                        y_sum += 1
                    elif j[k] == 5:
                        z_sum += 1
                u.append(u_sum), w.append(w_sum), x.append(x_sum), y.append(
                    y_sum
                ), z.append(z_sum)
    
            datas = {"N": n_ind_spec2, "Diagnostic": u, "Diag. or Partial": w, "Partial": x, "Part. or Uninformative": y, "Invalid": z}
            dc = pd.DataFrame(datas, index=spec_list2)   
            diagnostic.write(dc.to_string(index=False))
            return dc,spec_list2
    
    
    def plot_dc(self):
        list_var=self.list_all[1]
        cod_color=self.list_all[0]
        spec_list=self.spec_list2
        list_n_pos=self.list_n_pos
        """print HTML file with the resuls.              
        """     
        list_var = [x for x in list_var if x is not None]
        #    print list_var
        list_var2 = []
        list_n_pos = [
            list_n_pos[x] for x in range(len(cod_color)) if sum(cod_color[x]) != 0
        ]
        pos_cc = 0
        pos_cc_ls = []
        #    print cod_color
        for i in range(len(cod_color)):
            if sum(cod_color[i]) != 0:
                pos_cc_ls.append(pos_cc)
            pos_cc += 1
        #    print pos_cc_ls
        for i in range(len(list_var)):
            list_var2_tmp = []
            for j in pos_cc_ls:
                list_var2_tmp.append(list_var[i][j])
            list_var2.append(list_var2_tmp)
        list_var = list_var2
        #    print list_var
        cod_color = [
            cod_color[x] for x in range(len(cod_color)) if sum(cod_color[x]) != 0
        ]
        seq_snp = [[None] * len(list_var) for i in range(len(list_var[0]))]
        for i in range(len(list_var[0])):
            for j in range(len(list_var)):
                if cod_color[i][j] == 0:
                    seq_snp[i][j] = ""
                else:
                    seq_snp[i][j] = list_var[j][i][0]
        # print(seq_snp)
        # print(pos_cc_ls)
        # print(cod_color)
        # print(list_n_pos)
    
        #    Code for do the tables with the data prepared above
        n_div = 20 #number of columns in the table
        lists_data = [[] for i in range(int(ceil(len(pos_cc_ls) / float(n_div))))]
        div = int(ceil(len(pos_cc_ls) / float(n_div))) #number of tables necessaries if considering n_div columns
        n_col = 6 + ((len(spec_list) + 1) * div) #number of rows considering also the legend
        fct = 1 / float(n_col) #size of each row
        sz = round(1 - (6 * fct), 2)
        for i in range(int(ceil(len(pos_cc_ls) / float(n_div)))):  # 0 1 2
            if i + 1 < len(pos_cc_ls) / float(n_div):  # i < 2.1
                seq_snp_tmp = seq_snp[n_div * i : (n_div * (i + 1))]
            else:
                seq_snp_tmp = seq_snp[n_div * i :]
            if i + 1 < len(pos_cc_ls) / float(n_div):
                list_n_pos_tmp = list_n_pos[n_div * i : (n_div * (i + 1))]
            else:
                list_n_pos_tmp = list_n_pos[n_div * i :]
            if i + 1 < len(pos_cc_ls) / float(n_div):
                cod_color_tmp = cod_color[n_div * i : (n_div * (i + 1))]
            else:
                cod_color_tmp = cod_color[n_div * i :]
    
            cod_ini = []
            for j in range(len(cod_color_tmp[0])):
                cod_ini.append(0)
            cod_color2 = [cod_ini] + cod_color_tmp
            cod_color2 = [
                [
                    "red"
                    if x == 1
                    else "#ff5b00"
                    if x == 2
                    else "orange"
                    if x == 3
                    else "yellow"
                    if x == 4
                    else "#feffca"
                    if x == 4
                    else "white"
                    for x in l
                ]
                for l in cod_color2
            ]
            dict_color = {}
            dict_color["color"] = cod_color2
            head = ["<b>Group Name</b>"]
            data = {"Group_Name": spec_list}
            for k in range(len(seq_snp_tmp)):
                data[str(k)] = seq_snp_tmp[k]
            for l in list_n_pos_tmp:
                head.append("<b>" + str(l) + "</b>")
            df = pd.DataFrame(data)  #data frame with diagnostic characters
    
            val = []
            col = []
            for column in df:
                if column != "Group_Name":
                    col.append(column)
            col = [int(x) for x in col]
            col.sort()
            for m in col:
                val.append(df[str(m)])
            val = [df.Group_Name] + val
    
            width_1 = 200
            width_2 = 40
            yvalue = sz - ((len(spec_list) + 1) * fct)
            # print(yvalue)
            # print(sz)
            if yvalue < 0:
                yvalue = 0
            trace = go.Table(
                type="table",
                columnwidth=[width_1, width_2],
                domain=dict(
                    x=[
                        0,
                        float(width_1 + (width_2 * len(list_n_pos_tmp)))
                        / (width_1 + (width_2 * n_div)),
                    ],
                    y=[yvalue, sz - 0.02],
                ),
                header=dict(
                    values=head,
                    line=dict(color="#506784"),
                    fill=dict(color="white"),
                    align=["left", "center"],
                    font=dict(color="black", size=11),
                ),
                cells=dict(
                    values=val,
                    line=dict(color="#506784"),
                    fill=dict_color,
                    align=["left", "center"],
                    font=dict(color="black", size=11),
                ),
            )
    
            lists_data[i] = trace
    
            sz = sz - ((len(spec_list) + 1) * fct)
    
        data_lgnd = {
            "color": ["", "", "", "", ""],
            "Legend": [
                "Diagnostic Character",
                "Diagnostic or Partial Character",
                "Partial Character",
                "Partial or Uninformative Character",
                "Invalid Character",
            ],
        }
        df2 = pd.DataFrame(data_lgnd) #dataframe with legend
    
        trace_legend = go.Table(
            type="table",
            columnwidth=[width_2, 400],
            domain=dict(x=[0, 0.35], y=[round(1 - (6 * fct), 2), 1]),
            header=dict(
                values=["", "<b>Legend</b>"],
                line=dict(color="#506784"),
                fill=dict(color="white"),
                align=["left", "center"],
                font=dict(color="black", size=11),
            ),
            cells=dict(
                values=[df2.color, df2.Legend],
                line=dict(color="#506784"),
                fill=dict(
                    color=[
                        ["red", "#ff5b00", "orange", "yellow", "#feffca"],
                        "white",
                    ]
                ),
                align=["center"],
                font=dict(color="black", size=11),
            ),
        )
    
        height_per = (len(spec_list) + 1) * 30 * div + 6 * 60 # change for separate elements
        lists_data = [trace_legend] + lists_data
    
        layout2 = dict(
            width=950,
            height=height_per,
            autosize=False,
            title="Diagnostic Character detailed",
            margin=dict(t=100),
            showlegend=False,
        )
    
        fig = dict(data=lists_data, layout=layout2)
        plotly.offline.plot(
            fig,
            filename=os.path.join(self.path, "Diagnostic_character.html"),
            auto_open=False,
            config={"displayModeBar": False},
            show_link=False,
        )
        return (fig)

    def plot_legend(self):
        data_lgnd = {
            "color": ["", "", "", "", ""],
            "Legend": [
                "Diagnostic Character",
                "Diagnostic or Partial Character",
                "Partial Character",
                "Partial or Uninformative Character",
                "Invalid Character",
            ],
        }
        df2 = pd.DataFrame(data_lgnd) #dataframe with legend
        
        trace_legend = go.Table(
            type="table",
            columnwidth=[40, 400],
            domain=dict(x=[0, 0.35], y=[0, 1]),
            header=dict(
                values=["", "<b>Legend</b>"],
                line=dict(color="#506784"),
                fill=dict(color="white"),
                align=["left", "center"],
                font=dict(color="black", size=11),
            ),
            cells=dict(
                values=[df2.color, df2.Legend],
                line=dict(color="#506784"),
                fill=dict(
                    color=[
                        ["red", "#ff5b00", "orange", "yellow", "#feffca"],
                        "white",
                    ]
                ),
                align=["center"],
                font=dict(color="black", size=11),
            ),
        )
        
        layout2 = dict(
            height=250,
            autosize=True,
            margin=dict(l=10, r=0, t=20, b=0),
            showlegend=False,
                )
                
        fig = go.Figure(data=trace_legend, layout=layout2)
        fig.show()

    def plot_dcs(self):
        list_var=self.list_all[1]
        cod_color=self.list_all[0]
        spec_list=self.spec_list2
        list_n_pos=self.list_n_pos
        """print HTML file with the resuls.              
        """     
        list_var = [x for x in list_var if x is not None]
        #    print list_var
        list_var2 = []
        list_n_pos = [
            list_n_pos[x] for x in range(len(cod_color)) if sum(cod_color[x]) != 0
        ]
        pos_cc = 0
        pos_cc_ls = []
        #    print cod_color
        for i in range(len(cod_color)):
            if sum(cod_color[i]) != 0:
                pos_cc_ls.append(pos_cc)
            pos_cc += 1
        #    print pos_cc_ls
        for i in range(len(list_var)):
            list_var2_tmp = []
            for j in pos_cc_ls:
                list_var2_tmp.append(list_var[i][j])
            list_var2.append(list_var2_tmp)
        list_var = list_var2
        #    print list_var
        cod_color = [
            cod_color[x] for x in range(len(cod_color)) if sum(cod_color[x]) != 0
        ]
        seq_snp = [[None] * len(list_var) for i in range(len(list_var[0]))]
        for i in range(len(list_var[0])):
            for j in range(len(list_var)):
                if cod_color[i][j] == 0:
                    seq_snp[i][j] = ""
                else:
                    seq_snp[i][j] = list_var[j][i][0]
        # print(seq_snp)
        # print(pos_cc_ls)
        # print(cod_color)
        # print(list_n_pos)
    
        #    Code for do the tables with the data prepared above
        n_div = 20 #number of columns in the table
        lists_data = [[] for i in range(int(ceil(len(pos_cc_ls) / float(n_div))))]
        div = int(ceil(len(pos_cc_ls) / float(n_div))) #number of tables necessaries if considering n_div columns
        n_col = 6 + ((len(spec_list) + 1) * div) #number of rows considering also the legend
        fct = 1 / float(n_col) #size of each row
        sz = round(1 - (6 * fct), 2)
        for i in range(int(ceil(len(pos_cc_ls) / float(n_div)))):  # 0 1 2
            if i + 1 < len(pos_cc_ls) / float(n_div):  # i < 2.1
                seq_snp_tmp = seq_snp[n_div * i : (n_div * (i + 1))]
            else:
                seq_snp_tmp = seq_snp[n_div * i :]
            if i + 1 < len(pos_cc_ls) / float(n_div):
                list_n_pos_tmp = list_n_pos[n_div * i : (n_div * (i + 1))]
            else:
                list_n_pos_tmp = list_n_pos[n_div * i :]
            if i + 1 < len(pos_cc_ls) / float(n_div):
                cod_color_tmp = cod_color[n_div * i : (n_div * (i + 1))]
            else:
                cod_color_tmp = cod_color[n_div * i :]
    
            cod_ini = []
            for j in range(len(cod_color_tmp[0])):
                cod_ini.append(0)
            cod_color2 = [cod_ini] + cod_color_tmp
            cod_color2 = [
                [
                    "red"
                    if x == 1
                    else "#ff5b00"
                    if x == 2
                    else "orange"
                    if x == 3
                    else "yellow"
                    if x == 4
                    else "#feffca"
                    if x == 4
                    else "white"
                    for x in l
                ]
                for l in cod_color2
            ]
            dict_color = {}
            dict_color["color"] = cod_color2
            head = ["<b>Group Name</b>"]
            data = {"Group_Name": spec_list}
            for k in range(len(seq_snp_tmp)):
                data[str(k)] = seq_snp_tmp[k]
            for l in list_n_pos_tmp:
                head.append("<b>" + str(l) + "</b>")
            df = pd.DataFrame(data)  #data frame with diagnostic characters
    
            val = []
            col = []
            for column in df:
                if column != "Group_Name":
                    col.append(column)
            col = [int(x) for x in col]
            col.sort()
            for m in col:
                val.append(df[str(m)])
            val = [df.Group_Name] + val
    
            width_1 = 200
            width_2 = 40
            yvalue = sz - ((len(spec_list) + 1) * fct)
            # print(yvalue)
            # print(sz)
            if yvalue < 0:
                yvalue = 0
            trace = go.Table(
                type="table",
                columnwidth=[width_1, width_2],
                domain=dict(
                    x=[
                        0,
                        float(width_1 + (width_2 * len(list_n_pos_tmp)))
                        / (width_1 + (width_2 * n_div)),
                    ],
                    y=[yvalue, sz - 0.02],
                ),
                header=dict(
                    values=head,
                    line=dict(color="#506784"),
                    fill=dict(color="white"),
                    align=["left", "center"],
                    font=dict(color="black", size=11),
                ),
                cells=dict(
                    values=val,
                    line=dict(color="#506784"),
                    fill=dict_color,
                    align=["left", "center"],
                    font=dict(color="black", size=11),
                ),
            )
    
            lists_data[i] = trace
    
            sz = sz - ((len(spec_list) + 1) * fct)
    
    
        height_per = (len(spec_list) + 2) * 35 * div # change for separate elements
    
        layout2 = dict(
            width=950,
            height=height_per,
            autosize=False,
            # title="Diagnostic Character detailed",
            margin=dict(t=10),
            showlegend=False,
        )
    
        fig = go.Figure(data=lists_data, layout=layout2)
        # dict(data=lists_data, layout=layout2)
        fig.show()
    
def main(path, fasta, n_ind=3, gen=1, sp=2):
    """calculate and print diagnostic characters.
    Parameters
    ----------
    path: str
        The path to folder with the input file.
    fasta: str
        The name of fasta file.
    n_ind: int
        minimum number of individuals for species.        
    gen: int
        genus position on sequence names splited by "_".
    sp: int
        genus position on sequence names splited by "_".
    """
    # list_n_pos = snper(path, fasta)
    # n_spec = specier(os.path.join(path, "tmp/"), "snps_all.fasta", gen - 1, sp - 1)
    # spec_list = n_spec[1]
    # n_ind_spec = n_spec[2]
    # n_spec = n_spec[0]
    # list_all = diagnoser(path, n_spec, n_ind)
    # logging.info(("\n" + "### Summary of diagnostic character result ###\n"))
    datadc = diagnostic_character(path, fasta, n_ind, gen, sp)
    # datas = summary(path, list_all[0], spec_list, n_ind_spec, n_ind)
    # plot_dc(path, list_all[1], list_all[0], datas[1], list_n_pos)

    shutil.rmtree(os.path.join(path, "tmp"))  # will delete a directory and all its contents
    # os.remove() #will remove a file.
    # os.rmdir(path+'tmp') #will remove an empty directory.
    return datadc


if __name__ == "__main__":
    main()
#

# diagnoser("C:\Users\ramir\OneDrive\Python\Barcode example\Schizodon/",'SC_COI_ALINHADO.fasta',10,3)
# w.main()
