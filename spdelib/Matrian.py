# -*- coding: utf-8 -*-
"""
Created on Sun Apr 01 21:44:07 2018

@author: jolob
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#import numpy as np
from math import log, sqrt
import logging

class matrian: #lower matrix
    def __init__(self,path,fasta,gen,sp,distance):
        self.path=path
        self.fasta=open(path+fasta, newline='')
        self.lst_tot=self.lista_total()
        self.ncol=self.lst_tot[0][1]       
        self.Lname=self.name_sp(gen,sp)
        self.Lsp=self.listar_sp(sp)
        self.Lgen=self.listar_sp(gen)
        self.data=self.matrix(distance)

              
    def fasta_sp(self): #genera lista con nombres para hacer las otras listas
        self.fasta.seek(0)
        aln=self.fasta
#        print aln
        name=[]
        for line in aln:
            line = line.strip()
            if line.startswith(">"):
                name.append(line[1:])
#        print name        
        return name

    def fasta_seq(self): #genera lista de tuplas de nombre y sequencia
        self.fasta.seek(0)
        aln=self.fasta
        name=[]
        seq=[]
        D=[]
        for line in aln:
            line = line.strip()
            if line.startswith(">"):
                name.append(line[1:])
                D+=seq
                seq=['']
            elif not line:
                    pass
            elif line[0] != ">":
                    seq[0]+=(line.rstrip())
        D+=seq #apendea la ultima secuencia
        A=list(zip(name,D))                    
        return A
    
    def K2Pdistance(self,seq1,seq2):  #adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
        pairs = []
        for x in zip(seq1,seq2):
            if '-' not in x and 'N' not in x: pairs.append(x)     
        ts_count=0
        tv_count=0
        length = len(pairs)        
        transitions = [ "AG", "GA", "CT", "TC"]
        transversions = [ "AC", "CA", "AT", "TA","GC", "CG", "GT", "TG" ]    
        for (x,y) in pairs:
            if x+y in transitions: ts_count += 1 
            elif x+y in transversions: tv_count += 1
        p = float(ts_count) / length
        q = float(tv_count) / length
        try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )*100
        except ValueError: 
            logging.info("Tried to take log of a negative number")
            return None
        return d 



    def pdistance(self,seq1, seq2):  #adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py 
        p = 0    
        pairs = []    
        for x in zip(seq1,seq2):
            if '-' not in x and 'N' not in x: pairs.append(x)       
        for (x,y) in pairs:
            if x != y:   
                p += 1 
        length = len(pairs)   
        return (float(p) / length)*100

    
    def matrix(self,distance): #genera la matriz de distancia
        if distance == 'k':
            logging.info('Using k2p distance\n')
        else:
            logging.info('Using p-distance\n')
        m=[[None] * self.ncol for i in range(self.ncol)]
        A=self.fasta_seq()
        for i in range(self.ncol):
            for j in range(self.ncol):
#                if i==j or i>j:    ###Verificar si funciona!! y mejorar el tiempo
#                    continue
                if distance == 'k':
                    num=self.K2Pdistance(A[i][1],A[j][1])
                    if num==-0.0:
                        num=0.0
                    m[i][j]=(round(num,5)) #code by neguinha!                        
                if distance == 'p':
                    num=self.pdistance(A[i][1],A[j][1])
                    if num==-0.0:
                        num=0.0
                    m[i][j]=(round(num,5)) #code by neguinha!
        return m
                                
    def padron(self,i): #genera lista con padron '_'
        Padron=[0]
        Padro_pos=0
        for letra in i:
            Padro_pos+=1
            if letra=='_':
                Padron.append(Padro_pos)
        Padron.append(self.ncol)
        return Padron
    
    def listar_sp(self,a): #genera lista de especies o generos
        self.fasta.seek(0)
        L=[]
        sp=self.fasta_sp()
#        print sp
        for i in sp:           
            Padron=self.padron(i)
            ini=(Padron[a-1]+1)
            fin=(Padron[a]-1)           
            L.append(i[ini-1:fin])
        i1=''
        pos=0
        L1=[]
        for i in (L):
            pos+=1
            if i!=i1:
                L1.append(pos)
                i1=L[pos-1]
        L2=[]
        for i in range(len(L1)-1):
            L3=[L1[i],L1[i+1]-1]
            L2.append(L3)
        L2.append([L1[len(L1)-1],len(L)])        
        return L2      
        
    def lista_total(self):
        self.fasta.seek(0)
        T=self.fasta_sp()
        L=[[1,len(T)]]
        return L
        
    def name_sp(self,a,b): #genera lista de nombres de especies
        self.fasta.seek(0)
        L=[]
        sp=self.fasta_sp()
        for i in sp:
            Padron=self.padron(i)
            ini=(Padron[a-1]+1)
            fin2=(Padron[b]-1)
            L.append(i[ini-1:fin2])  
        Lname=[]
        [Lname.append(key) for key in L if key not in Lname]
        return Lname
        
    def conc_intra(self,ini,fin): #concatena intra dentro de una sp
        if fin>ini:
            x=ini-1
            y=ini-1
            z=ini-1
            conc=[]
            while x<fin and z<fin:
                while y<=z and y!=x:
                    conc.append(self.data[x][y])
                    y+=1
                z+=1
                y=ini-1    
                x+=1
            return conc
        else:
            return None
            
    def conc_inter(self,ini,fin,ini2,fin2): #conc inter de um sp por genero
        if ini==ini2 and fin==fin2:
            return None
        x=ini-1
        y=ini2-1
        conc=[]
        while x<fin:
            while y<ini-1:
                conc.append(self.data[x][y])
                y+=1
            y=ini2-1   
            x+=1
        x=fin
        y=ini-1
        while x<fin2:
            while y<fin:
                conc.append(self.data[x][y])           
                y+=1
            y=ini-1   
            x+=1   
        return conc   
                             
    def intra(self,ini,fin): #da la media intra para un dado grupo
        L=self.conc_intra(ini,fin)
        if L==None:
            return None
        v=0
        for i in L:
            v+=i
        m=(float(v))/(len(L))
        return m
               
    def intra_all(self): #da media intra normalizada de varios grupos #used
        m=0
        n=0
        for i in self.Lsp:
            if i[1]>i[0]:
                m+=self.intra(i[0],i[1])
                n+=1
        M=m/n
        return M
              
    def inter(self,ini,fin,ini2,fin2): #da media interespecifica por genero
        if ini==ini2 and fin==fin2:
            return None
        L=self.conc_inter(ini,fin,ini2,fin2)
        v=0
        for i in L:
            v+=i
        m=(float(v))/(len(L))
        return m
        
    def inter_all(self): #da media intergrupo normalizada #used
        m=0
        n=0
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    if self.inter(e[0],e[1],g[0],g[1])!=None:
                        m+=self.inter(e[0],e[1],g[0],g[1])
                        n+=1
        M=m/n
        return M
        
    def interALL_all(self): #da media inter general normalizada #used
        m=0
        n=0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    if self.inter(e[0],e[1],g[0],g[1])!=None:
                        m+=self.inter(e[0],e[1],g[0],g[1])
                        n+=1
        M=m/n
        return M

    def inter_intra_all_name(self): #da media inter e intra de todos los grupos y nombre #used
        II=[]
        x=0
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    name=self.Lname[x]
                    x+=1
                    inter=self.inter(e[0],e[1],g[0],g[1])
                    intra=self.intra(e[0],e[1])
                    intera=[name,inter,intra]
                    II.append(intera)
        return II

    def interALL_intra_all_name(self): #da media inter general e intra de todos los grupos y nombre #used
        II=[]
        x=0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    name=self.Lname[x]
                    x+=1
                    inter=self.inter(e[0],e[1],g[0],g[1])
                    intra=self.intra(e[0],e[1])
                    intera=[name,inter,intra]
                    II.append(intera)
        return II
        
                        
    def max_intra(self,ini,fin): #da maxima dentro de un grupo
        L=self.conc_intra(ini,fin)
        if L==None:
            return None
        v=max(L)
        return v
        
    def min_intra(self,ini,fin): #da minima dentro de un grupo
        L=self.conc_intra(ini,fin)
        if L==None:
            return None
        v=min(L)
        return v
        
    def min_inter(self,ini,fin,ini2,fin2): #da minima inter de um grupo
        L=self.conc_inter(ini,fin,ini2,fin2)
        if L==None:
            return None
        v=min(L)
        return v
        
    def max_inter(self,ini,fin,ini2,fin2): #da maxima inter de um grupo
        L=self.conc_inter(ini,fin,ini2,fin2)
        if L==None:
            return None
        v=max(L)
        return v
        
    def min_inter_pos(self,ini,fin,ini2,fin2): #posicion en lista conc_inter del minima inter de um grupo #used
        L=self.conc_inter(ini,fin,ini2,fin2)
        if L==None:
            return None
        vpos=L.index(min(L))
        return vpos        
    
    def min_all2(self): #da min intergrupos #used
        MM=[]
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    mi=self.min_inter(e[0],e[1],g[0],g[1])
                    if mi!=None:
                        MM.append(mi)
        return MM
        
    def max_inter_all(self): #da max intergrupos #used
        MM=[]
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    mi=self.max_inter(e[0],e[1],g[0],g[1])
                    MM.append(mi)
        return MM

    def max_inter_ALL(self): #da max inter sin importar los grupos #used
        MM=[]
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    mi=self.max_inter(e[0],e[1],g[0],g[1])
                    MM.append(mi)
        return MM
        
    def min_ALL(self): #da min inter general sin importar los grupos #used
        MM=[]
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    mi=self.min_inter(e[0],e[1],g[0],g[1])
                    MM.append(mi)
        return MM
        
    def max_all(self): #da max intra de todos los grupos #used
        MM=[]
        for i in self.Lsp:
            ma=self.max_intra(i[0],i[1])
            MM.append(ma)
        return MM
        
    def min_intra_all(self): #da min intra de todos los grupos #used
        MM=[]
        for i in self.Lsp:
            ma=self.min_intra(i[0],i[1])
            if ma!=None:
                MM.append(ma)
        return MM
        
    def min_max_ALL_name(self): #da min general y max intra con nombre de sp       #used
        MM=[]
        x=0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    no=self.Lname[x]
                    x+=1
                    mi=self.min_inter(e[0],e[1],g[0],g[1])
                    ma=self.max_intra(e[0],e[1])
                    mima=[no,mi,ma]
                    MM.append(mima)
        return MM
        
    def min_ALL_pos(self): #da min general y max con nombre de sp      #used 
        MMpos=[]
        MM_num_indv_x_esp=[]
        MM_pos_name_indv=[]
        MM_name_indv_NN=[]
        i=0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    pos_mi=self.min_inter_pos(e[0],e[1],g[0],g[1])
                    MMpos.append(pos_mi)   
                num_indv_x_esp=self.Lsp[i][1]-self.Lsp[i][0]+1
                MM_num_indv_x_esp.append(num_indv_x_esp)
                if MMpos[i] < (MM_num_indv_x_esp[i]*(self.Lsp[i][0]-1)):
                    pos_name_indv=(MMpos[i]%(self.Lsp[i][0]-1))
                else:
                    pos_name_indv=(MMpos[i]//MM_num_indv_x_esp[i])+(MM_num_indv_x_esp[i])    
                MM_pos_name_indv.append(pos_name_indv) 
                i+=1
        for j in MM_pos_name_indv:   
            name_indv_min_inter=self.fasta_seq()[j][0]  
            MM_name_indv_NN.append(name_indv_min_inter)   
        return MM_name_indv_NN
        
    def min_media_max_intra(self): #used
        new_max = [x if x!=None else 0 for x in self.max_all()] 
        MX=max(new_max)
        new_min = [x if x!=None else 0 for x in self.min_intra_all()]         
        MI=min(new_min)
        med=self.intra_all()
        return MI,med,MX   
        
    def min_media_max_inter(self): #used
        new_max = [x if x!=None else 0 for x in self.max_inter_all()] 
        MX=max(new_max)
        new_min = [x if x!=None else 0 for x in self.min_all2()]  
        MI=min(new_min)
        med=self.inter_all()
        return MI,med,MX              

    def min_media_max_tot(self): #used
        MX=max(self.max_inter_ALL())
        MI=min(self.min_ALL())
        med=self.interALL_all()
        return MI,med,MX 
            
    def cont_intra_all(self): #cont freq intra de todos los grupos
        cont=[]
        for i in self.Lsp:
            if i[1]>i[0]:
                cont+=(self.conc_intra(i[0],i[1]))
        cont.sort()
        #cont2={x:cont.count(x) for x in set(cont)}
        #return cont2
        return cont    
        
    def cont_inter_all(self): #cont freq inter de todos los grupos dentro de su genero
        cont_inter=[]
        for g in self.Lgen: #Lgenlst_tot
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    if (self.conc_inter(e[0],e[1],g[0],g[1]))!=None:
                        cont_inter+=(self.conc_inter(e[0],e[1],g[0],g[1]))
        cont_inter.sort()
        #cont2={x:(cont_inter.count(x)/2) for x in set(cont_inter)}
        #return cont2
        return cont_inter
        
    def inter_geo(self): #da freq unico inter de todos los grupos dentro de su grupo mayor
        cont_inter=[]
        x=0
        for g in self.Lgen: #Lgenlst_tot
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    if (self.conc_inter(e[0],e[1],g[0],g[1]))!=None:
                        inter=(self.conc_inter(e[0],e[1],g[0],g[1]))
                        inter.sort()
                        inter2=[]
                        for e in inter:
                            if e not in inter2:
                                inter2.append(e)
                        no=self.Lname[x]
                        x+=1
                        inno=[no,inter2] 
                        cont_inter.append(inno)           
        return cont_inter

    def cont_inter_geo(self): #cont freq unico inter de todos los grupos dentro de su grupo mayor
        cont_inter=[]
        for g in self.Lgen: #Lgenlst_tot
            for e in self.Lsp:
                if e[1]<=g[1] and e[0]>=g[0]:
                    if (self.conc_inter(e[0],e[1],g[0],g[1]))!=None:
                        inter=(self.conc_inter(e[0],e[1],g[0],g[1]))
                        inter.sort()
                        inter2=[]
                        for e in inter:
                            if e not in inter2:
                                inter2.append(e)
                        cont_inter.append(inter2)
        y=0
        cont2=[]
        cont3=[]
        for e in range(len(self.Lsp)):
            if y==e:
                cont2.append(cont_inter[y])
                y+=2
        for e in cont2:
            for i in e:
                cont3.append(i)                                    
        cont3.sort()
        cont4={x:(cont3.count(x)) for x in set(cont3)}
        return cont4
  

    def plot_freq(self): #Barcoding gap graph 
        ter=self.cont_inter_all()
        tra=self.cont_intra_all()
        newBins_tra=len(set(tra))//3
        newBins_ter=len(set(ter))//3
        f, ax = plt.subplots(3, 1, sharex='col', sharey='all')
        sns.distplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance',ax=ax[1])
        sns.distplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance',ax=ax[0])
        sns.distplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance',ax=ax[2])
        sns.distplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance',ax=ax[2])   
        ax[0].set_title('DNA Barcoding gap')
        ax[0].set_ylabel('# of taxon pairs')
        ax[1].set_ylabel('# of taxon pairs')
        ax[2].set_ylabel('# of taxon pairs')
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        ax[2].set_xlabel('Genetic Distance')
        f.savefig(self.path+'barcoding_gap.pdf')
#        plt.show(f)
#        plt.close(f)
        plt.clf()




        
    def med_ind_sp(self): #media de individuso por sp
        m=0
        t=0
        for i in self.Lsp:
            m+=i[1]-i[0]+1
            t+=1
        return float(m/t)

    def analyze(self): #hace analisis en bloque
        logging.info('Summary table (Name, mean intra, max intra, NN, distance to NN) in percentage')   
        a,b,c,d,e=[],[],[],[],[]  
        II=self.inter_intra_all_name()
        mima=self.min_max_ALL_name()
        name_NNindv=self.min_ALL_pos()        
        for i in range(len(II)):
            a.append(II[i][0])
            b.append(II[i][2])
            c.append(mima[i][2])
            d.append(name_NNindv[i])
            e.append(mima[i][1])
        summ=pd.DataFrame()
        summ['Name']=a
        summ['Mean']=b
        summ['Max']=c
        summ['NN']=d
        summ['DtoNN']=e
        logging.info(summ.to_string(index=False))          
        logging.info('')
        logging.info('min interspecific and max intraspecific by group')    
        mima=self.min_max_ALL_name()
        x=[]
        y=[]
        y2=[]
        z=[]       
        for i in mima:
            x.append(i[1])
            y2.append(i[2])
            if i[2] == None:
                y.append(0.0)
            else:
                y.append(i[2])
            z.append(i[0])
        df=pd.DataFrame()
        df['name']=z
        df['inter']=x
        df['intra']=y2
        df['intra2']=y
        df1=df[['name','inter','intra']].copy()
        logging.info(df1.to_string(index=False))
        
        logging.info('') 
        tra=[]
        ter=[]
        total=[]
        title=['minimum', 'mean', 'maximum']
        L=self.min_media_max_intra()
        for i in L:
            tra.append(i)      
        M=self.min_media_max_inter()
        for i in M:
            ter.append(i)
        tot=self.min_media_max_tot()
        for i in tot:
            total.append(i)
        tab=pd.DataFrame()
        tab['']=title
        tab['intra']=tra
        tab['inter']=ter
#        tab['total']=total
        logging.info(tab.transpose().to_string(header=False))



        sns.lmplot('intra2', 'inter', 
           data=df, 
           fit_reg=False) 
        plt.title('Maximum intraspecific vs Minimum to NN')
        plt.xlabel('Maximum intraspecific')
        plt.ylabel('Minimum to NN')
        z=[max(x),max(y)]
        plt.axis([0,max(z)+1,0,max(z)+1])
        lims= [0,max(z)+1]
        plt.plot(lims, lims, ':k')        
        plt.savefig(self.path+'min_max.pdf')
#        plt.show()
        plt.clf()
 
        
        self.plot_freq()
        
    def geoanalyze(self): #hace analisis geobarcode
        logging.info('valores inter')
        ter=self.inter_geo()
        for i in ter:        
            logging.info(i)
        ter2=self.cont_inter_geo()
        for k,v in list(ter2.items()):
            logging.info(k,v)
            
def main(path,fasta,gen,sp,distance,out_name=None,n=False):
#    targets = logging.StreamHandler(sys.stdout), logging.FileHandler(path+'SPdel_output.log') # from https://stackoverflow.com/questions/24204898/python-output-on-both-console-and-file/24206109#24206109
#    logging.basicConfig(format='%(message)s', level=logging.INFO, handlers=targets)
    tmp=matrian(path,fasta,gen,sp,distance)
    tmp.analyze()
    if n==True:
        List=tmp.name_sp(gen,sp)
        with open(out_name,'w+') as nominal_file:
            for i in List:
                nominal_file.write(i+'\n')   

    
    
    
    
if __name__ == "__main__":
    main()
       