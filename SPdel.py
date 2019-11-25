#SPdel v 1.0
try:
    import sys
    import os
    from Bio import SeqIO
    from ete3 import TextFace, RectFace, Tree, TreeStyle, NodeStyle
    from ete3.treeview.faces import add_face_to_node
    import Diagnoser
    import Matrian
    import random
    import colorsys
except ImportError:
	print("Please install all dependences first.")
	sys.exit()

def check_list(path,fasta,fileList):
    with open(path+fasta, newline='') as fasta, open(path+fileList, newline='') as filelist:
        n_line=0
        n_seq=0
        for line in filelist:
            n_line+=1
        for line in fasta:
            if line.startswith('>'):
                n_seq+=1
        if n_seq==n_line:
            return True
        
def rename_NomList(path,fasta,NomList,tree=None):
    if not os.path.exists(path+'tmp_file/'):
        os.makedirs(path+'tmp_file/')
    newpath=path+'tmp_file/'
    with open(path+fasta, newline='') as fasta, open(path+NomList, newline='') as newnames, open(newpath+'fasta_renamed.fasta', 'w+') as newfasta:
        if tree!=None:
            with open(path+tree, newline='') as tree, open(newpath+'tree_renamed.nwk', 'w+') as newtree:
                for line in tree:
                    new_line_tree=line
                for line in fasta:
                    if line.startswith('>'):
                        new_name= newnames.readline().strip().replace(' ','_')
                        newline=line.replace('>','').split(' ')
        #                print new_name
                        newfasta.write('>'+new_name+'_'+newline[0])
                        new_line_tree=new_line_tree.replace(newline[0].strip(),(new_name+'_'+newline[0].strip()))
                    else:
        #                print line
                        newfasta.write(line)
                newtree.write(new_line_tree)
        else:
                for line in fasta:
                    if line.startswith('>'):
                        new_name= newnames.readline().strip().replace(' ','_')
                        newline=line.replace('>','').split(' ')
        #                print new_name
                        newfasta.write('>'+new_name+'_'+newline[0])
                    else:
        #                print line
                        newfasta.write(line)
                        
def rename_MOTUList(path,fasta,List,NomList,outfile):  
    if not os.path.exists(path+'tmp_file/'):
        os.makedirs(path+'tmp_file/')
    newpath=path+'tmp_file/'                      
    with open(path+fasta, newline='') as fasta, open(path+NomList, newline='') as newnames, open(path+List, newline='') as PL, open(newpath+outfile, 'w+') as newPL:
        for line in PL:
            new_line_PL=line
#            print(new_line_PL)
            fasta.seek(0)
            newnames.seek(0)
            for line in fasta:
                if line.startswith('>'):
                    new_name= newnames.readline().strip().replace(' ','_')
                    newline=line.replace('>','').split(' ')
                    new_line_PL=new_line_PL.replace(newline[0].strip(),(new_name+'_'+newline[0].strip()))
#                print(new_line_PL)
            newPL.write(new_line_PL)                       
                        
def random_color(h=None, l=None, s=None):
    """ returns the RGB code of a random color. Hue (h), Lightness (l)
    and Saturation (s) of the generated color could be fixed using the
    pertinent function argument.  """
    def rgb2hex(rgb):
        return '#%02x%02x%02x' % rgb
    def hls2hex(h, l, s):
        return rgb2hex( tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(h, l, s))))

    if not h:
        h = random.random()
    if not s: 
        s = 0.5
    if not l:
        l = 0.5
    return hls2hex(h, l, s)

class sorting_fasta:
    def __init__(self,path,fasta,typeCODE='VER',outname='Nominal/sorted.fasta'):    
        self.path=path
        self.fasta=open(path+fasta, newline='')
        self.outname=outname
        self.typeCODE=typeCODE

    def check_fasta(self):
        print('checking if sequences are aligned')
        self.fasta.seek(0)
        handle = self.fasta
        n=len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
#        print n
        for seq in SeqIO.parse(handle, "fasta"):
            if(len(seq))==n:
#                print "OK"
                continue
            else:
#                print 'error'
                return False
#        print 'ok'
        
    def describe_fasta(self):
        self.fasta.seek(0)
        handle = self.fasta
        n=len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
        n_seq=0
        for seq in SeqIO.parse(handle, "fasta"):
                n_seq+=1
        print('Fasta file with',n_seq,'sequences and',n,'base pairs')
        

    def stop_codon(self):    #modified from https://stackoverflow.com/questions/34009041/python-code-to-find-coding-dna-with-start-and-stop-codons      
        print('checking stop codons using '+self.typeCODE+' genetic code' )
#        print(self.typeCODE)
        if self.typeCODE=='VER':
            codon_list = ["AGA", "TAA", "TAG", "AGG"] #The Vertebrate Mitochondrial Code (transl_table=2)
            codon_list_RC = ["TCT", "TTA", "CTA", "CCT"]
        elif self.typeCODE=='INV':
            codon_list = ["TAA", "TAG"] #The Invertebrate Mitochondrial Code (transl_table=5)
            codon_list_RC = ["TTA", "CTA"] 
        self.fasta.seek(0)
        handle = self.fasta
        k=[0,1,2]
        H=[-1,-2,-3]
        for seq in SeqIO.parse(handle, "fasta"):
            codonPass=[]
            for i in k:
                n=i
#                print(n)
                found_codon_positions = []
                lenseq=len(seq)
                while n < lenseq-2:
                    possible_codon = seq[n:n+3]
                    if possible_codon.seq in codon_list:
                        found_codon_positions.append(n+1)
                    n += 3
                if found_codon_positions!=[]:
#                    print('found codons at indices {}'.format(found_codon_positions)+' at '+seq.id)                 
                    codonPass.append(1)
                else:
                    codonPass.append(0)
            for i in H:
                n=i
#                print(n)
                found_codon_positions = []
                lenseq=len(seq)
                while n*-1 < lenseq-2:
                    if n==0:
                        possible_codon = seq[n-3:]
                    else:
                        possible_codon = seq[n-3:n]                    
                    if possible_codon.seq in codon_list_RC:
                        found_codon_positions.append(n+1)
                    n += -3
                if found_codon_positions!=[]:
#                    print('found codons at indices {}'.format(found_codon_positions)+' at '+seq.id)                 
                    codonPass.append(1)
                else:
                    codonPass.append(0)                    
#            print(codonPass)
            if 0 in codonPass:
                pass
            else:
                print('WARNING: Any ORF without stop codon were found in sequence '+seq.id)
           
        
            
    def sort_fasta(self):
        print('sorting fasta')
        self.fasta.seek(0)
        handle = self.fasta
        l = SeqIO.parse(handle, "fasta")
        sortedList = [f for f in sorted(l, key=lambda x : x.id)]
        tmp=self.path+self.outname
        file_sorted= open(tmp, 'w+')
        for s in sortedList:
#            print s.description
            file_sorted.write('>'+s.description+'\n')
#            return s.description
#            print str(s.seq)
            file_sorted.write(str(s.seq)+'\n')
        print('checking if sorting was ok')
        handle.seek(0)
        file_sorted.seek(0)
        n_aln=0
        n_aln_sorted=0
        for seq in SeqIO.parse(handle, "fasta"):
            n_aln+=1
#        print n_aln
        for seq in SeqIO.parse(file_sorted, "fasta"):
            n_aln_sorted+=1
#        print n_aln_sorted
        if n_aln == n_aln_sorted:
            print('Fasta file sorted successfuly')
            pass
        else:
            self.sort_fasta(self)
#        handle.close()
        file_sorted.close()
        return file_sorted
    
    def check_pattern(self):
        self.fasta.seek(0)
        for line in self.fasta:
            if line.startswith('>'):
                pattern=line.split('_')
                if len(pattern)<3:
                    return False
     
class check_tree:
    def __init__(self,path,fasta,tree):    
        self.path=path
        self.fasta=open(path+fasta, newline='')
        self.tree=Tree(path+tree)
        
    def check(self):
        self.fasta.seek(0)
        t=self.tree
        name_leaf=[]
        for leaf in t:
            name_leaf.append(leaf.name)
#        print len(name_leaf)
        name_fasta=[]
        for seq in SeqIO.parse(self.fasta, "fasta"):
            name_fasta.append(seq.id)
#        print name_leaf
#        print name_fasta
#        print len(name_leaf)
#        print len(name_fasta)
        if len(name_leaf) != len(name_fasta):
            print('Different number of individuals in tree and fasta')
            print('Names found in tree: '+'%s' % ', '.join(map(str, list(set(name_leaf) - set(name_fasta)))))
            print('Names found in fasta: '+'%s' % ', '.join(map(str, list(set(name_fasta) - set(name_leaf)))))              
            return False
        comp=set(name_leaf) & set(name_fasta)
#        print len(comp)
        if len(comp) != len(name_leaf):
            print('Names found in tree: '+'%s' % ', '.join(map(str, list(set(name_leaf) - set(name_fasta)))))
            print('Names found in fasta: '+'%s' % ', '.join(map(str, list(set(name_fasta) - set(name_leaf)))))              
            return False
      
    def root_midle(self):
        t = self.tree
#        print t
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
        return t
        
    
class MOTU_PTP:
    def __init__(self,path,fasta,outfile,outptp='PTP'):
        if not os.path.exists(path+'PTP/'):
            os.makedirs(path+'PTP/')
        self.path=path
        self.fasta_sorted=open(path+fasta, newline='')
        self.outfile=outfile
        self.outptp=outptp
        self.fasta_PTP=self.MOTU_renameFasta()
              
    def MOTU_listPTP(self):
        listPTP=open(self.path+self.outfile)
        listPTP.seek(0)
        name=[]
        for line in listPTP:
            line = line.strip()
            if not line.strip(): continue
            if line.startswith("#"):
                pass
            elif line.startswith("Species"):
                pass
            else:
                name.append(line)
        return name
    
    def print_MOTU_PTP(self):
        listPTP=open(self.path+self.outfile)
        listPTP.seek(0)
        listPTP_txt=open(self.path+'PTP/'+self.outptp+'.MOTU_list.txt', 'w+')
        for line in listPTP:
            line = line.strip()
            if not line.strip(): continue
            if line.startswith("#"):
                pass
            elif line.startswith("Species"):
                print('')
                line = line.replace("Species","MOTU")
                print(line)
                listPTP_txt.write(line+'\n')
            else:
                line = line.replace(",",", ")
                print(line)
                listPTP_txt.write(line+'\n')

    def MOTU_renameFasta(self):
        self.fasta_sorted.seek(0)
        newfasta= open(self.path+'PTP/'+'PTP_MOTU.fasta', 'w+')
        a=0
        for i in self.MOTU_listPTP():   
            sps= i
            sp=sps.split(',')
            sp.sort()
#            print(sp)
            for j in sp:
#                print j
                for line in self.fasta_sorted:
#                    if line.startswith('>') and str(j)+'\n' in line:
                    if '>'+str(j) in line:                        
                        new_name= '>MOTU_'+str(a+1)+'_'+str(j)+'\n'
                        newfasta.write(new_name)
                        while True:
                            seq=next(self.fasta_sorted,'')
                            if seq=='':
                                break
                            elif seq[0]!='>':
                                newfasta.write(seq)
                            else:
                                break
                    else:
                        pass
                self.fasta_sorted.seek(0)               
            a+=1
#checking fasta size number of sequences
        self.fasta_sorted.seek(0)
        newfasta.seek(0)
        n_aln=0
        n_aln_sorted=0
        for seq in SeqIO.parse(self.fasta_sorted, "fasta"):
            n_aln+=1
#        print(n_aln)
        for seq in SeqIO.parse(newfasta, "fasta"):
            n_aln_sorted+=1
#        print(n_aln_sorted)
        if n_aln == n_aln_sorted:
#            print('Fasta file sorted successfuly')
            pass
        else:
            print('error renaming fasta')
#cheking fasta for bp size
        newfasta_check=sorting_fasta(self.path+'PTP/','PTP_MOTU.fasta')   #Checking fasta PROVISIONALLY!!!
        if newfasta_check.check_fasta() is False:
            sys.exit("Error: different length sequences, please check your alignment")
        #
        newfasta.seek(0)
        return newfasta    
    
class MOTU_bPTP:
    def __init__(self,path,fasta,outfile,outbptp='bPTP'):
        if not os.path.exists(path+'bPTP/'):
            os.makedirs(path+'bPTP/')
        self.path=path
        self.fasta_sorted=open(path+fasta, newline='')
        self.outptp=outbptp
        self.outfile=outfile
        self.fasta_bPTP=self.MOTU_renameFasta2()
     
    def MOTU_listbPTP(self):
        listPTP=open(self.path+self.outfile)
        listPTP.seek(0)
        name=[]
        for line in listPTP:
            line = line.strip()
            if not line.strip(): continue
            if line.startswith("#"):
                pass
            elif line.startswith("Species"):
                pass
            else:
                name.append(line)
        return name

    def print_MOTU_bPTP(self):
        listPTP=open(self.path+self.outfile)
        listPTP.seek(0)
        listPTP_txt=open(self.path+'bPTP/'+self.outptp+'.MOTU_list.txt', 'w+')
        for line in listPTP:
            line = line.strip()
            if not line.strip(): continue
            if line.startswith("#"):
                pass
            elif line.startswith("Species"):
                print('')
                line = line.replace("Species","MOTU")
                print(line)
                listPTP_txt.write(line+'\n')
            else:
                line = line.replace(",",", ")
                print(line)
                listPTP_txt.write(line+'\n')

    def MOTU_renameFasta2(self):
        self.fasta_sorted.seek(0)
        newfasta= open(self.path+'bPTP/'+'bPTP_MOTU.fasta', 'w+')
        a=0
        for i in self.MOTU_listbPTP():   
            sps= i
            sp=sps.split(',')
            sp.sort()
    #        print sp
            for j in sp:
                for line in self.fasta_sorted:
                    if '>'+str(j) in line:
                        new_name= '>MOTU_'+str(a+1)+'_'+str(j)+'\n'
    #                    print new_name
                        newfasta.write(new_name)
                        while True:
                            seq=next(self.fasta_sorted,'')
                            if seq=='':
                                break
                            elif seq[0]!='>':
                                newfasta.write(seq)
                            else:
                                break
                    else:
                        pass
                self.fasta_sorted.seek(0)               
            a+=1
#        self.fasta_sorted.close()
#        newfasta.close()
#checking fasta size
        self.fasta_sorted.seek(0)
        newfasta.seek(0)
        n_aln=0
        n_aln_sorted=0
        for seq in SeqIO.parse(self.fasta_sorted, "fasta"):
            n_aln+=1
#        print n_aln
        for seq in SeqIO.parse(newfasta, "fasta"):
            n_aln_sorted+=1
#        print n_aln_sorted
        if n_aln == n_aln_sorted:
#            print 'Fasta file sorted successfuly'
            pass
        else:
            print('error renaming fasta')
#            self.MOTU_renameFasta2()
#
        newfasta.seek(0)
        return newfasta    
    
class MOTU_GMYC:
    def __init__(self,path,fasta,outfile):
        if not os.path.exists(path+'GMYC/'):
            os.makedirs(path+'GMYC/')
        self.path=path
        self.fasta_sorted=open(path+fasta, newline='')
        self.outfile=outfile
        self.fasta_gmyc=self.MOTU_renameFasta_gmyc()
        
    def MOTU_listGMYC(self):
        listGMYC=open(self.path+self.outfile)
        listGMYC.seek(0)
        name=[]
        for line in listGMYC:
            line = line.strip()
            line = line.replace(" ","")
            line = line[:-1]
            if not line.strip(): continue
            if line.startswith("Species"):
                pass
            else:
                name.append(line)
#        print name
        return name
    
    def print_MOTU_GMYC(self):
        listGMYC=open(self.path+self.outfile)
        listGMYC.seek(0)
        listGMYC_txt=open(self.path+'GMYC/'+'MOTU_list.txt', 'w+')
        for line in listGMYC:
            line = line.strip()
            if not line.strip(): continue
            if line.startswith("Species"):
                print('')
                line = line.replace("Species","MOTU")
                print(line)
                listGMYC_txt.write(line+'\n')
            else:
                line = line[:-1]
                print(line)
                listGMYC_txt.write(line+'\n')

    def MOTU_renameFasta_gmyc(self):
        self.fasta_sorted.seek(0)
        newfasta= open(self.path+'GMYC/'+'GMYC_MOTU.fasta', 'w+')
        a=0
        for i in self.MOTU_listGMYC():   
            sps= i
            sp=sps.split(',')
            sp.sort()
#            print sp
            for j in sp:
                for line in self.fasta_sorted:
                    if '>'+str(j) in line:
                        new_name= '>MOTU_'+str(a+1)+'_'+str(j)+'\n'
    #                    print new_name
                        newfasta.write(new_name)
                        while True:
                            seq=next(self.fasta_sorted,'')
                            if seq=='':
                                break
                            elif seq[0]!='>':
                                newfasta.write(seq)
                            else:
                                break
                    else:
                        pass
                self.fasta_sorted.seek(0)               
            a+=1
#        self.fasta_sorted.close()
#        newfasta.close()
#checking fasta size
        self.fasta_sorted.seek(0)
        newfasta.seek(0)
        n_aln=0
        n_aln_sorted=0
        for seq in SeqIO.parse(self.fasta_sorted, "fasta"):
            n_aln+=1
#        print n_aln
        for seq in SeqIO.parse(newfasta, "fasta"):
            n_aln_sorted+=1
#        print n_aln_sorted
        if n_aln == n_aln_sorted:
#            print 'Fasta file sorted successfuly'
            pass
        else:
            print('error renaming fasta')
#            self.MOTU_renameFasta_gmyc(self)
#
        newfasta.seek(0)
        return newfasta  

class MOTU_BIN:
    def __init__(self,path,fasta,BinList=None,):
        if not os.path.exists(path+'BIN/'):
            os.makedirs(path+'BIN/')
        self.path=path+'BIN/'
        self.fasta=open(path+fasta, newline='')
        self.BinList=path+BinList
        self.fasta_bin=self.MOTU_renameFasta_bin()

    def print_MOTU_bin(self):
        with open(self.path+'BIN_MOTU_sorted.fasta') as listbin, open(self.path+'MOTU_list.txt', 'w+') as listBIN_txt:
#        listbin=open(self.path+'BIN_MOTU_sorted.fasta')
#        listBIN_txt=open(self.path+'MOTU_list.txt', 'w+')
            listbin.seek(0)
            name_list=[]
            motu_list=[]
#            a=0
#            for line in listbin:
#                a+=1
#            print a
            listbin.seek(0)
            for line in listbin:
                if line.startswith('>'):
                    line.strip()
                    bins=line.split('_')
#                    print bins ###########
#                    print bins[1] ###########
                    motu_name=bins[1]
                    if motu_name not in motu_list:
                        motu_list.append(motu_name)
            listbin.seek(0)    
            for line in listbin:
                if line.startswith('>'):
                    line.strip()
#                    print line
                    line=line[:-1]
                    x=line.split('_')  ###gran descubrimiento!
#                    print x
                    name_list.append(x)
            for i in motu_list:
                member=[]
                print(('MOTU '+i))
                listBIN_txt.write('MOTU '+i+'\n')
                for j in name_list:
#                    print j
                    if i==j[1]:  
                        a=2
                        y=len(j)
                        ind=''
                        while a<y:
                            ind=ind+j[a]+'_'
                            a+=1
                        member.append(ind[:-1])
                print('%s' % ', '.join(map(str, member)))
                listBIN_txt.write('%s' % ', '.join(map(str, member))+'\n')
                print('')

    def MOTU_renameFasta_bin(self):
        self.fasta.seek(0)
        with open(self.BinList) as newnames, open(self.path+'BIN_MOTU.fasta', 'w+') as newfasta:
            for line in self.fasta:
                if line.startswith('>'):
                    bin_name=newnames.readline().strip()
                    new_name= '>MOTU_'+bin_name+'_'+line[1:]#+'\n'
    #                print new_name
                    newfasta.write(new_name)
                else:
                    newfasta.write(line)
        newfasta2=sorting_fasta(self.path,'BIN_MOTU.fasta',outname='BIN_MOTU_sorted.fasta')   #Checking fasta
        if newfasta2.check_fasta() is False:
            sys.exit("Error: different length sequences, please check your alignment")            
        newfasta2.sort_fasta()          #sorting fasta
        return newfasta2



class MOTU_X:
    def __init__(self,path,fasta,folder,BinList=None):
        if not os.path.exists(path+folder):
            os.makedirs(path+folder)
        self.path=path+folder+'/'
        self.fasta=open(path+fasta, newline='')
        self.BinList=path+BinList
        self.fasta_X=self.MOTU_renameFasta_X()

    def print_MOTU_X(self):
        with open(self.path+'X_MOTU_sorted.fasta') as listbin, open(self.path+'MOTU_list.txt', 'w+') as listBIN_txt:

            listbin.seek(0)
            name_list=[]
            motu_list=[]

            listbin.seek(0)
            for line in listbin:
                if line.startswith('>'):
                    line.strip()
                    bins=line.split('_')
                    motu_name=bins[1]
                    if motu_name not in motu_list:
                        motu_list.append(motu_name)
            listbin.seek(0)    
            for line in listbin:
                if line.startswith('>'):
                    line.strip()
#                    print line
                    line=line[:-1]
                    x=line.split('_') 
#                    print x
                    name_list.append(x)
            for i in motu_list:
                member=[]
                print(('MOTU '+i))
                listBIN_txt.write('MOTU '+i+'\n')
                for j in name_list:
#                    print j
                    if i==j[1]:  
                        a=2
                        y=len(j)
                        ind=''
                        while a<y:
                            ind=ind+j[a]+'_'
                            a+=1
                        member.append(ind[:-1])
                print('%s' % ', '.join(map(str, member)))
                listBIN_txt.write('%s' % ', '.join(map(str, member))+'\n')
                print('')

    def MOTU_renameFasta_X(self):
        self.fasta.seek(0)
        with open(self.BinList) as newnames, open(self.path+'X_MOTU.fasta', 'w+') as newfasta:
            for line in self.fasta:
                if line.startswith('>'):
                    bin_name=newnames.readline().strip()
                    new_name= '>MOTU_'+bin_name+'_'+line[1:]#+'\n'
    #                print new_name
                    newfasta.write(new_name)
                else:
                    newfasta.write(line)
        newfasta2=sorting_fasta(self.path,'X_MOTU.fasta',outname='X_MOTU_sorted.fasta')   #Checking fasta
        if newfasta2.check_fasta() is False:
            sys.exit("Error: different length sequences, please check your alignment")            
        newfasta2.sort_fasta()          #sorting fasta
        return newfasta2
    
    
    
class MOTU_nominal:
    def __init__(self,path,fasta,nominalList=None,gen=1,sp=2):
        self.path=path
        self.fasta=open(path+fasta, newline='')
        if nominalList is None:
            self.Lname=self.print_MOTU_nonminal(gen,sp)
        else:
            self.nominalList=path+nominalList  #PARA OPCION DE LISTA CON NOMINALES!!!
        
        
    def print_MOTU_nonminal(self,a,b):
        with open(self.path+'tmp_file/nominal_list.txt', newline='') as listnominal_txt, open(self.path+'Nominal/nominal_list_sp.txt', 'w+') as allList:
            self.fasta.seek(0)
            name=[]
            n=0
            for line in self.fasta:
                line = line.strip()
                if line.startswith(">"):
                    name.append(line[1:])
            for i in listnominal_txt:
                i=i.strip()
                n+=1
#                print('Species '+str(n)+' '+i)
                allList.write('MOTU '+str(n)+' '+i+'\n')
                member=[]
                i=i+'_'
                for j in name:
                    if i in j:
                        member.append(j)
#                print ('%s' % ', '.join(map(str, member))+'\n')
                allList.write('%s' % ', '.join(map(str, member))+'\n')                        
                                
class Compare:
    def __init__(self,path,fasta,CompList):
        if not os.path.exists(path+'Compare/'):
            os.makedirs(path+'Compare/')
        self.path=path
        self.CompList=CompList.split(',')
        self.compared=self.compare()
        self.summary=self.summary_motu()
        self.fasta=open(path+fasta, newline='')
#        self.fasta_comp=self.MOTU_renameFasta_Compare()
#        self.comp_warn=self.Compare_warning()

        
    def compare(self):
        with open(self.path+'tmp_file/ALL_MOTUs.txt','w+') as compare_file:
            lists=[]
            motus_name=[]
            
            if 'n' in self.CompList and os.path.exists(self.path+'Nominal/nominal_list_sp.txt'):
                listNOM=open(self.path+'Nominal/nominal_list_sp.txt')
                for line in listNOM:
                    if line.startswith('MOTU'):
                        newline=line.strip().split()[:2]
                        newline.append('(NOM)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(';'))
            if 'p' in self.CompList and os.path.exists(self.path+'PTP/'+'PTP'+'.MOTU_list.txt'):
                listPTP=open(self.path+'PTP/'+'PTP'+'.MOTU_list.txt')
                for line in listPTP:
                    if line.startswith('MOTU'):
                        newline=line.strip().split()[:2]
                        newline.append('(PTP)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(';'))
            if 't' in self.CompList and os.path.exists(self.path+'bPTP/'+'bPTP'+'.MOTU_list.txt'):
                listbPTP=open(self.path+'bPTP/'+'bPTP'+'.MOTU_list.txt')
                for line in listbPTP:
                    if line.startswith('MOTU'): 
                        newline=line.strip().split()[:2]
                        newline.append('(bPTP)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(";"))
            if 'g' in self.CompList and os.path.exists(self.path+'GMYC/'+'MOTU_list.txt'):
                listGMYC=open(self.path+'GMYC/'+'MOTU_list.txt')
                for line in listGMYC:
                    if line.startswith('MOTU'): 
                        newline=line.strip().split()[:2]
                        newline.append('(GMYC)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(";"))
            if 'b' in self.CompList and os.path.exists(self.path+'BIN/'+'MOTU_list.txt'):
                listBIN=open(self.path+'BIN/'+'MOTU_list.txt')
                for line in listBIN:
                    if line.startswith('MOTU'): 
                        newline=line.strip().split()[:2]
                        newline.append('(BIN)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(";"))
            if 's' in self.CompList and os.path.exists(self.path+'OT/'+'OT_MOTU_list.txt'):
                listGMYC=open(self.path+'OT/'+'OT_MOTU_list.txt')
                for line in listGMYC:
                    if line.startswith('MOTU'): 
                        newline=line.strip().split()[:2]
                        newline.append('(OT)')
                        motus_name.append(newline)
                    else:
                        line = line.replace(", ",";")
                        lists.append(line.strip().split(";"))
                        
#            if self.XList!=None:                                    
#                XList=self.XList.split(',')
            for i in self.CompList:
                if i != 'n' and 't' and 'p' and 'g' and 's' and 'b':
                    i_folder=i.split('.')[0]                        
                    if os.path.exists(self.path+i_folder+'/'+'MOTU_list.txt'):
                        listX=open(self.path+i_folder+'/'+'MOTU_list.txt')
                        for line in listX:
                            if line.startswith('MOTU'): 
                                newline=line.strip().split()[:2]
                                newline.append('('+i_folder+')')
                                motus_name.append(newline)
                            else:
                                line = line.replace(", ",";")
                                lists.append(line.strip().split(";"))

                                      
            concordants=[]
            used=[]
            pos_names=[]
            for i in range(len(lists)):
                elem_sort=lists[i]
                elem_sort.sort()
                if elem_sort in used:
                    continue
                pos_names_for=[]
                used.append(elem_sort)
                pos_names_for.append(i)
                concordant=[]
                concordant.append(lists[i])
                for j in range(len(lists)):
                    if i==j or i>j:
                        continue
                    h=lists[i]
                    k=lists[j]
                    if len(set(h))==len(set(k)) and len((set(h).intersection(set(k))))==len(set(h)):
                        concordant.append(k)
                        pos_names_for.append(j)
                concordants.append(concordant)
                pos_names.append(pos_names_for)
    #        a=1
            for i in range(len(concordants)):
    #            print ('')                    
    #            print 'Concordant MOTU',a
    #            a+=1
                motus_con=[]
                for j in pos_names[i]:
                    con='%s' % ' '.join(map(str, motus_name[j]))
                    motus_con.append(con)
    #            print ('%s' % ' & '.join(map(str, motus_con)))
                compare_file.write(('%s' % ' & '.join(map(str, motus_con)))+'\n')
    #            print ('%s' % ', '.join(map(str, concordants[i][0])))
                compare_file.write(('%s' % ', '.join(map(str, concordants[i][0])))+'\n')
    
    
    def summary_motu(self):
        n_analysis=0
        if 'p' in self.CompList and os.path.exists(self.path+'PTP/'+'PTP'+'.MOTU_list.txt'):
            n_analysis+=1
        if 't' in self.CompList and os.path.exists(self.path+'bPTP/'+'bPTP'+'.MOTU_list.txt'):
            n_analysis+=1
        if 'g' in self.CompList and os.path.exists(self.path+'GMYC/'+'MOTU_list.txt'):
            n_analysis+=1
        if 'b' in self.CompList and os.path.exists(self.path+'BIN/'+'MOTU_list.txt'):
            n_analysis+=1
        if 's' in self.CompList and os.path.exists(self.path+'OT/'+'OT_MOTU_list.txt'):
            n_analysis+=1
        for i in self.CompList:
            if i != 'n' and 't' and 'p' and 'g' and 's' and 'b':
#            for i in XList:    
                i_folder=i.split('.')[0]                        
                if os.path.exists(self.path+i_folder+'/'+'MOTU_list.txt'):
                    n_analysis+=1
            
    #    print n_analysis
        if os.path.exists(self.path+'tmp_file/ALL_MOTUs.txt'):
            with open(self.path+'tmp_file/ALL_MOTUs.txt', newline='') as compare_file, open(self.path+'Compare/Summary_MOTUs.txt', 'w+') as summary_file, open(self.path+'Compare/Compare_MOTU_list.txt', 'w+') as out_list:
                num_motu=1
                lists_TC=[]
                motus_name_TC=[]
                lists_MC=[]
                motus_name_MC=[]
                lists_TD=[]
                motus_name_TD=[]
                lists_MD=[]
                motus_name_MD=[]
                
                a='########################################\n'
                if 'n' in self.CompList and os.path.exists(self.path+'Nominal/nominal_list_sp.txt'):
                    for line in compare_file:
                        if line.startswith('MOTU') and '(NOM)' in line:
                            n_used=0
                            n_used=len(line.strip().split("&"))-1
#                            print(n_used,n_analysis,n_analysis/float(2)) ###############
                            if n_used==n_analysis:
                                motus_name_TC.append(line)
                                lists_TC.append(next(compare_file))
                            if n_used<n_analysis and n_used>=n_analysis/float(2):
                                motus_name_MC.append(line)
                                lists_MC.append(next(compare_file))                            
                        if line.startswith('MOTU') and not '(NOM)' in line:
                            n_used=0
                            n_used=len(line.strip().split("&"))
#                            print(n_used,n_analysis,n_analysis/float(2)) ###############
                            if n_used==n_analysis:
                                motus_name_TD.append(line)
                                lists_TD.append(next(compare_file))
                            if n_used<n_analysis and n_used>=n_analysis/float(2):
                                motus_name_MD.append(line)
    #                            print line
                                lists_MD.append(next(compare_file))                                             
                    
                    print(('\n'+'### MOTUs totally concordant with taxonomy ###\n'))
                    summary_file.write('\n'+a+'MOTUs totally concordant with taxonomy\n'+a+'\n')
                    for i in range(len(motus_name_TC)):
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_TC[i].strip()+']')
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_TC[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print(lists_TC[i])
                        summary_file.write(lists_TC[i]+'\n')
                        out_list.write(lists_TC[i])
                        num_motu+=1
                    print(('\n'+'### MOTUs mostly concordant with taxonomy ###\n'))
                    summary_file.write('\n'+a+'MOTUs Mostly concordant with taxonomy\n'+a+'\n')
                    for i in range(len(motus_name_MC)):
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_MC[i].strip()+']')
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_MC[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print(lists_MC[i])
                        summary_file.write(lists_MC[i]+'\n')
                        out_list.write(lists_MC[i])
                        num_motu+=1                            
                    print(('\n'+'### MOTUs totally discordant with taxonomy ###\n'))
                    summary_file.write('\n'+a+'MOTUs totally discordant with taxonomy\n'+a+'\n')
                    for i in range(len(motus_name_TD)):
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_TD[i].strip()+']')
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_TD[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print(lists_TD[i])
                        summary_file.write(lists_TD[i]+'\n')
                        out_list.write(lists_TD[i])
                        num_motu+=1
                    print(('\n'+'### MOTUs mostly discordant with taxonomy ###\n'))
                    summary_file.write('\n'+a+'MOTUs mostly discordant with taxonomy\n'+a+'\n')
                    for i in range(len(motus_name_MD)):
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_MD[i].strip()+']')
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_MD[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print(lists_MD[i])
                        summary_file.write(lists_MD[i]+'\n')
                        out_list.write(lists_MD[i])
                        num_motu+=1
                       
                else:
                    for line in compare_file:
                        if line.startswith('MOTU'):
                            n_used=0
                            n_used=len(line.strip().split("&"))
                            print(n_used,n_analysis,n_analysis/float(2)) ###############
                            if n_used==n_analysis:
                                motus_name_TC.append(line)
                                lists_TC.append(next(compare_file))
                            if n_used<n_analysis and n_used>=n_analysis/float(2):
                                motus_name_MC.append(line)
                                lists_MC.append(next(compare_file))
 
                    print(('\n'+'### MOTUs totally concordant ###\n'))
                    summary_file.write('\n'+a+'MOTUs totally concordant\n'+a+'\n')
                    for i in range(len(motus_name_TC)):
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_TC[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_TC[i].strip()+']')
                        summary_file.write(lists_TC[i]+'\n')
                        out_list.write(lists_TC[i])
                        print(lists_TC[i])
                        num_motu+=1
                    print(('\n'+'### MOTUs mostly concordant ###\n'))
                    summary_file.write('\n'+a+'MOTUs mostly concordant\n'+a+'\n')
                    for i in range(len(motus_name_MC)):
                        print('Consensus MOTU '+str(num_motu)+' ['+motus_name_MC[i].strip()+']')
                        summary_file.write('Consensus MOTU '+str(num_motu)+' ['+motus_name_MC[i].strip()+']\n')
                        out_list.write('MOTU '+str(num_motu)+'\n')
                        print(lists_MC[i])
                        summary_file.write(lists_MC[i]+'\n')
                        out_list.write(lists_MC[i])
                        num_motu+=1

    def Compare_warning(self):
        lists_w_temp=[]
        lists_W=[]
        motus_name_W=[]
        motus_memb_W=[]
        a='########################################\n'
        with open(self.path+'Compare/Summary_MOTUs.txt', 'r+') as summary_file:
            summary_file.seek(0)
            for line in summary_file:                        
                if line.startswith('Consensus'):
                    motus_name_W.append(line.strip())
                    member=next(summary_file)
                    lists_w_temp.append(member.strip().split(', '))
                    motus_memb_W.append(member)
            for i in range(len(lists_w_temp)):   #MC vs MC
                for j in range(len(lists_w_temp)):
                    if i!=j and i>j:
                        if set(lists_w_temp[i]) & set(lists_w_temp[j]) != set([]):
                            lists_W.append([i,j]) 
                            
            if lists_W != []:
                print(('\n'+'### WARNING: overlapping MOTUS  ###\n'))
                summary_file.write('\n'+a+'WARNING: overlapping MOTUS\n'+a+'\n')
                n=0
                for i in lists_W:
                    n+=1
                    print(('\n'+'Overlapping consensus MOTUs '+str(n)+'\n'))
                    summary_file.write('\n'+'Overlapping consensus MOTUs '+str(n)+'\n')
                    print(motus_name_W[i[0]].strip())
                    summary_file.write(motus_name_W[i[0]])
                    print(motus_memb_W[i[0]])
                    summary_file.write(motus_memb_W[i[0]]+'\n')
                    print(motus_name_W[i[1]].strip())
                    summary_file.write(motus_name_W[i[1]])
                    print(motus_memb_W[i[1]])
                    summary_file.write(motus_memb_W[i[1]]+'\n')

                warn_chos=str(input('Do you want to choice between the overlapping MOTUs before calculate DNA barcoding statistics? (Y/N)').strip())
                while warn_chos!='Y' and warn_chos!='y'and warn_chos!='n'and warn_chos!='N':
                    print('Error: Answer must be Y or N')
                    warn_chos=str(input('Do you want to chose between the overlapping MOTUs before calculate DNA barcoding statistics? (Y/N)').strip())
                if warn_chos == 'Y' or warn_chos == 'y':
                    rmvd,acpt,num_del=[],[],[]
                    for i in lists_W:
                        motuA=motus_name_W[i[0]].strip()
                        motuB=motus_name_W[i[1]].strip()
                        text_cho='1)'+motuA+' or 2)'+motuB+' (1/2)'
                        motu_cho=int(input(text_cho))
                        while motu_cho!=1 and motu_cho!=2:
                            print('Error: Answer must be 1 or 2')
                            motu_cho=int(input(text_cho))
                        print('\nMOTU chosen '+motus_name_W[i[motu_cho-1]]+'\n')
                        summary_file.write('MOTU chosen '+motus_name_W[i[motu_cho-1]]+'\n')
                        if motu_cho-1 == 0:
                            motu_del=1
                        else:
                            motu_del=0
                        rmvd.append(motus_name_W[i[motu_del]])
                        acpt.append(motus_name_W[i[motu_cho-1]])
                        num_del.append(motu_del)
                        
                    while set(rmvd) & set(acpt) != set([]):
                        print('You cant remove and acepted the same MOTU!!!!\n')
                        
                        rmvd,acpt=[],[]
                        for i in lists_W:
                            motuA=motus_name_W[i[0]].strip()
                            motuB=motus_name_W[i[1]].strip()
                            text_cho='1)'+motuA+' or 2)'+motuB+' (1/2)'
                            motu_cho=int(input(text_cho))
                            while motu_cho!=1 and motu_cho!=2:
                                print('Error: Answer must be 1 or 2')
                                motu_cho=int(input(text_cho))
                            print('\nMOTU chosen '+motus_name_W[i[motu_cho-1]]+'\n')
                            summary_file.write('MOTU chosen '+motus_name_W[i[motu_cho-1]]+'\n')
                            if motu_cho-1 == 0:
                                motu_del=1
                            else:
                                motu_del=0
                            rmvd.append(motus_name_W[i[motu_del]])
                            acpt.append(motus_name_W[i[motu_cho-1]])
                            
    
                    for ind,i in enumerate(lists_W):
                        with open(self.path+'Compare/Compare_MOTU_list.txt', 'r+') as out_list:
                            list_out=out_list.readlines()
                        with open(self.path+'Compare/Compare_MOTU_list.txt', 'w+') as out_list:
                            for ind2,line in enumerate(list_out):
                                if line.startswith('MOTU'):
                                    new_line=list_out[ind2+1]
                                    if new_line == motus_memb_W[i[num_del[ind]]]:
                                        pass
                                    else:
                                        out_list.write(line)
                                        out_list.write(new_line)
    
                elif warn_chos == 'N' or warn_chos == 'n':
                    print('You need to choice a consensus MOTUS to calculate distances and graphs')
                    exit()


    def MOTU_renameFasta_Compare(self):
        self.fasta.seek(0)
        with open(self.path+'Compare/Compare_MOTU_list.txt', newline='') as Comp_file, open(self.path+'Compare/Compare_MOTU.fasta', 'w+') as newfasta:
            for line in Comp_file:
                if line.startswith('MOTU'):
                    new_name=line.replace(' ','_').strip()
                else:
                    line=line.replace(' ','').strip().split(',')
                    for i in line:
                        self.fasta.seek(0)
                        s=SeqIO.parse(self.fasta, "fasta")
                        for seq in s:
#                            print seq.id
                            if seq.id==i:
                                newfasta.write('>'+new_name+'_'+i+'\n')
#                                print seq.seq
                                newfasta.write(str(seq.seq)+'\n')
                        
                
class plot_compare_tree:  #cambiar nombres y cursiva, mostrar texto nominal?,
    def __init__(self,path,tree,names=True):    #names=False for hide names leaf
        self.path=path
        self.tree_loaded=Tree(path+tree)
        self.names=names
        self.dic_color=self.dic_motu_color()
        self.plot_tree_con()

    def dic_motu_color(self):
        color_used=['#e6194b','#3cb44b','#ffe119','#0082c8','#f58231','#911eb4','#46f0f0','#f032e6','#d2f53c','#fabebe','#008080','#e6beff','#aa6e28','#fffac8','#800000','#aaffc3','#808000','#ffd8b1','#000080','#808080','#FFFFFF','#000000']
        with open(self.path+'tmp_file/ALL_MOTUs_con.txt','w+') as motus_file:
            motus_con_file=open(self.path+'Compare/'+'Compare_MOTU_list.txt')    
            motus_file2=open(self.path+'tmp_file/ALL_MOTUs.txt')
            for line in motus_file2:
                if line.startswith('MOTU'):
                    motus=next(motus_file2)
                    motus2 = motus.replace(", ",";")
                    h=motus2.strip().split(';')
                    con_be=0
                    for j in motus_con_file:
                        if j.startswith('MOTU'):
                            cons=next(motus_con_file)
                            cons = cons.replace(", ",";")
                            k=cons.strip().split(';')                        
                            if len(set(h))==len(set(k)) and len((set(h).intersection(set(k))))==len(set(h)):
                                motus_file.write(line.strip()+' & (CONSENSUS)'+'\n')
                                motus_file.write(motus)
                                con_be=+1

                    if con_be==0:
                        motus_file.write(line)
                        motus_file.write(motus)                                
                    motus_con_file.seek(0)
                    
        with open(self.path+'tmp_file/ALL_MOTUs_con.txt') as motus_file:
            Ana_name=['NOM']
            for line in motus_file:
                if line.startswith('MOTU'):
                    motus=line.strip().split('&')
                    for i in range(len(motus)):
                        name=motus[i].replace(')','')
                        name=name.replace(' ','')    
                        name=name.split('(')
                        if name[1] not in Ana_name:
                            Ana_name.append(name[1])
            Ana_name.remove('CONSENSUS')
            Ana_name.append('CONSENSUS')

        with open(self.path+'tmp_file/ALL_MOTUs_con.txt') as motus_file:            
            n=0
            dicts=dict((name,{}) for name in Ana_name)
#            print(dicts)
            for line in motus_file:
                if line.startswith('MOTU'):
                    if n<22:
                        color=color_used[n]
                        n+=1
                    else:
                        color=random_color()
                        while True:
                            if color not in color_used:
                                color_used.append(color)
                                break
                            else:
                                color=random_color()                    
                    prevline=line
                    motus=next(motus_file)
                    motus = motus.replace(", ",";")
                    lists=motus.strip().split(';')

                    for name in Ana_name:
                        name2='('+name+')'
                        if name2 in prevline:
                            for i in lists:
                                dicts[name][i]=color
                                
        return dicts
                           
    def layout_tree_compare(self,node):
        node.collapsed = False
        n=0
        if node.is_leaf():
            name = node.get_leaf_names()[0]
            for k in self.dic_color:
                n+=1
                add_face_to_node(RectFace(10,10, 'FFFFFF', self.dic_color[k][name]), node, n, aligned=True)

            
    def plot_tree_con(self):          
        ts = TreeStyle()
        ts.margin_left=15
        ts.margin_bottom=15
        ts.margin_right=15
        ts.margin_right=15
#        ts.show_border=True
        ts.show_leaf_name=self.names
        ts.layout_fn = self.layout_tree_compare
        ts.title.add_face(TextFace("Species delimitation comparision", ftype='Arial', fsize=10), column=0)
        n=0
        for k in self.dic_color:
            F=TextFace(k, fsize=7)
            F.rotation=270            
            n+=1
            ts.aligned_header.add_face(F, column=n)            
        nstyle=NodeStyle()
        nstyle['size']=0
        nstyle['vt_line_width']=1
        nstyle['hz_line_width']=1
        for n in self.tree_loaded.traverse():
            n.set_style(nstyle)
#        self.tree_loaded.show(tree_style=ts)      #delete?? kernel died every twice 
        self.tree_loaded.render(self.path+"Compare/compare_tree.pdf", tree_style=ts,dpi=300)


def print_options():
    print('')
    print('SPdel v1.0 - Species delimitation and statistics for DNA Barcoding data sets')
    print('')
    print('The sequences name should be separate for "_" (e.g. Genus_species_individual) or use -N option for rename sequences')
    print('')        
    print('usage: ./SPdel.py path_to_files/ fasta_file -a n -P PTP_File -t tree_file -X MOTUList1.txt,MOTUList2.txt -C p,MOTUList1,SC_MOTUList2 -code VER')
    print('usage: ./SPdel.py path_to_files/ fasta_file -a n -distance p')
    print('usage: ./SPdel.py path_to_files/ fasta_file -a -P PTP_File -G GMYC_File -T bPTP_File -t tree_file -B BIN_file')
    print('usage: ./SPdel.py path_to_files/ fasta_file -a ptg\n')
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
    print("    -D           For Diagnostic character analysis.\n")
    print("    -C           Specify the type of analisys to be compared, include n for nominal, p for PTP, t for bPTP, b for BIN, and any filename used in X option for external MOTU list.\n")    
    print("    -code        Specify the genetic code used to test stop codon, VER or INV.\n")    
    print("Options for nominal: \n")
    print('    -gen         Position of the genus name in the sequence name when split by "_" (default=1).\n')
    print('    -sp          Position of the species name in the sequence name when split by "_" (default=2)\n')
    print("Options for diagnostic character: \n")
    print("    -n_ind       Minimum number of individuals for species to be considered in the diagnostic character analysis (default=3)\n")

    
def main(argu=None):
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
    path=sys.argv[1]
    fasta=sys.argv[2]
    a=''
    gen=1
    sp=2
    tree=None
    dis='k'
    BinList=None
    XList=None
    NomList=None
    CompList=None
    PTPList=None
    GMYCList=None
    bPTPList=None
    n_ind=3
    CODE=None
    
    for i in range(len(sys.argv)):
        if sys.argv[i] =='-n':
            a+='n'
        elif sys.argv[i] =='-t':
            i=i+1
            tree=sys.argv[i]
        elif sys.argv[i] =='-distance':
            i=i+1
            dis=sys.argv[i]
        elif sys.argv[i] =='-B':
            i=i+1
            a+='b'
            BinList=sys.argv[i]
        elif sys.argv[i] =='-D':
            a+='d'
        elif sys.argv[i] =='-P':
            i=i+1
            a+='p'
            PTPList=sys.argv[i]
        elif sys.argv[i] =='-G':
            i=i+1
            a+='g'
            GMYCList=sys.argv[i]
        elif sys.argv[i] =='-T':
            i=i+1
            a+='t'
            bPTPList=sys.argv[i]            
        elif sys.argv[i] =='-X':
            i=i+1
            a+='x'            
            XList=sys.argv[i]
        elif sys.argv[i] =='-N':
            i=i+1
            NomList=sys.argv[i]        
        elif sys.argv[i] =='-gen':
            i=i+1
            gen=sys.argv[i]
        elif sys.argv[i] =='-sp':
            i=i+1
            sp=sys.argv[i]            
        elif sys.argv[i] =='-n_ind':
            i=i+1
            n_ind=sys.argv[i]
        elif sys.argv[i] =='-code':
            i=i+1
            CODE=sys.argv[i]
        elif sys.argv[i] =='-C':
            i=i+1
            CompList=sys.argv[i] 
           
    tex='#####################\n'
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    if not os.path.exists(path+'tmp_file/'):
        os.makedirs(path+'tmp_file/')
    if not os.path.exists(path+'Nominal/'):
        os.makedirs(path+'Nominal/')
    if any(x in a for x in ('n','p','t','g','b','x')):
        if NomList!=None:
            check_file=check_list(path,fasta,NomList)
            if check_file==True:
                pass
            else:
                sys.exit('Different number of individuals in Nominal List provided')                
            print('Renaming sequences using nominal list provided')
            if PTPList!=None:
                rename_MOTUList(path,fasta,PTPList,NomList,'MOTU_PTP_List_renamed.txt')
                PTPList='tmp_file/MOTU_PTP_List_renamed.txt'            
            if bPTPList!=None:
                rename_MOTUList(path,fasta,bPTPList,NomList,'MOTU_bPTP_List_renamed.txt')
                bPTPList='tmp_file/MOTU_bPTP_List_renamed.txt' 
            if GMYCList!=None:
                rename_MOTUList(path,fasta,GMYCList,NomList,'MOTU_GMYC_List_renamed.txt')
                GMYCList='tmp_file/MOTU_GMYC_List_renamed.txt'                 
            if tree!=None:
                rename_NomList(path,fasta,NomList,tree)
                tree='tmp_file/tree_renamed.nwk'
                fasta='tmp_file/fasta_renamed.fasta'
            else:
                rename_NomList(path,fasta,NomList)
                fasta='tmp_file/fasta_renamed.fasta'                
            check=sorting_fasta(path,'tmp_file/fasta_renamed.fasta',typeCODE=CODE)   #Checking fasta
            if check.check_fasta() is False:
                sys.exit("Error: different length sequences, please check your alignment")            
            check.describe_fasta()  #describing fasta
            if CODE!=None:
                check.stop_codon()
            check.sort_fasta()          #sorting fasta
        else:
            check=sorting_fasta(path,fasta,typeCODE=CODE)   #Checking fasta
            if check.check_fasta() is False:
                sys.exit("Error: different length sequences, please check your alignment")            
            check.describe_fasta()  #describing fasta
            if CODE!=None:            
                check.stop_codon()
            check.sort_fasta()          #sorting fasta
    if 'n' in a:
        if check.check_pattern()==False:
            sys.exit('Please rename your file using "genera_especies_individual" format or provide a nominal list file')
        print(('\n'+tex+'Nominal species\n'+tex))
        Matrian.main(path+'Nominal/','sorted.fasta',gen,sp,dis,path+'tmp_file/nominal_list.txt',n=True)
        MOTU_nominal(path,'Nominal/sorted.fasta')
        if 'd' in a:
            Diagnoser.main(path+'Nominal/','sorted.fasta',n_ind)
    if tree!=None:
        ch_tr=check_tree(path,fasta,tree)
        if ch_tr.check() is False: 
           sys.exit("Error: different names in tree, please check your tree and alignment")
#            new_tree=ch_tr.root_midle() ####### test
    if 'p' in a:
        print(('\n'+tex+'PTP MOTUs\n'+tex))
#        print tree
        sys.argv=[sys.argv[0]]
        PTP_analize=MOTU_PTP(path,'Nominal/sorted.fasta',PTPList)   #tree=new_tree) test rooted tree
        print('\n### Delimited MOTUs ###')
        PTP_analize.print_MOTU_PTP()
        print('')
        Matrian.main(path+'PTP/','PTP_MOTU.fasta',gen,sp,dis)
        if 'd' in a:
            Diagnoser.main(path+'PTP/','PTP_MOTU.fasta',n_ind)
    if 't' in a:
        print(('\n'+tex+'bPTP MOTUs\n'+tex))         
        sys.argv=[sys.argv[0]]
        bPTP_analize=MOTU_bPTP(path,'Nominal/sorted.fasta',bPTPList)
        print('')
        bPTP_analize.print_MOTU_bPTP()
        print('\n### Delimited MOTUs ###')
        Matrian.main(path+'bPTP/','bPTP_MOTU.fasta',gen,sp,dis)
        if 'd' in a:
            Diagnoser.main(path+'bPTP/','bPTP_MOTU.fasta',n_ind)
    if 'g' in a:
        print(('\n'+tex+'GMYC MOTUs\n'+tex))
        sys.argv=[sys.argv[0]]
        gmyc_analize=MOTU_GMYC(path,'Nominal/sorted.fasta',GMYCList)
        print('\n### Delimited MOTUs ###')
        gmyc_analize.print_MOTU_GMYC()
        print('')
        Matrian.main(path+'GMYC/','GMYC_MOTU.fasta',gen,sp,dis)
        if 'd' in a:
            Diagnoser.main(path+'GMYC/','GMYC_MOTU.fasta',n_ind)
    if 'b' in a:
        if BinList==None:
            sys.exit('Please provide a BIN List')
        check_file=check_list(path,fasta,BinList)
        if check_file==True:
            pass
        else:
            sys.exit('Different number of individuals in BIN List provided') 
        print(('\n'+tex+'BIN MOTUs\n'+tex))
        sys.argv=[sys.argv[0]]
        bin_analize=MOTU_BIN(path,fasta,BinList=BinList)
        print('\n### Delimited MOTUs ###')
        bin_analize.print_MOTU_bin()
        print('')
        Matrian.main(path+'BIN/','BIN_MOTU_sorted.fasta',gen,sp,dis)
        if 'd' in a:
            Diagnoser.main(path+'BIN/','BIN_MOTU_sorted.fasta',n_ind)
    if 'x' in a:
        if XList==None:
            sys.exit('Please provide a MOTUs List')
        XList2=XList.split(',')
        for i in XList2:
            check_file=check_list(path,fasta,i)
            if check_file==True:
#                print(i)
                pass
            else:
                sys.exit('Different number of individuals in MOTUs List provided')

            i_folder=i.split('.')[0]        
            print(('\n'+tex+i_folder+' MOTUs\n'+tex))
            sys.argv=[sys.argv[0]]
            
            X_analize=MOTU_X(path,fasta,i_folder,BinList=i)
            print('\n### Delimited MOTUs ###')
            X_analize.print_MOTU_X()
            print('')
            Matrian.main(path+i_folder+'/','X_MOTU_sorted.fasta',gen,sp,dis)
            if 'd' in a:
                Diagnoser.main(path+i_folder+'/','X_MOTU_sorted.fasta',n_ind)

    if CompList!=None:     
        print(('\n'+tex+'Concordant MOTUs\n'+tex))
        comp_analize=Compare(path,fasta,CompList)
        comp_analize.Compare_warning()
        comp_analize.MOTU_renameFasta_Compare()
        if tree!=None:
            plot_compare_tree(path,tree)
        Matrian.main(path+'Compare/','Compare_MOTU.fasta',gen,sp,dis)
        if 'd' in a:
            Diagnoser.main(path+'Compare/','Compare_MOTU.fasta',n_ind)
        
        
#    sys.stdout = orig_stdout
#    f.close()
        
        
if __name__ == "__main__":
#    basepath = os.path.dirname(os.path.abspath(__file__))
    main()    

   
#optimizar matrian
    #revisar total distancia
    #usar split para patron

#tk

#output print to file #HTML?

#path fasta

# test N or gaps
#testar reenraizar a la mitad?

#separar opcion compare
    #casos en que mi sobre posicion no se da bien....cada analisis una delimit
    