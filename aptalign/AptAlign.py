import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptamat2.0'))
sys.path.append(root_path)
root_path = os.path.abspath(os.path.join(current_dir, '..','clustering'))
sys.path.append(root_path)
#import clustering_AptaMat as CLAM
import AptaMat2 as AF
import time
import argparse
import multiprocessing as mp
# import matplotlib.pyplot as plt
import numpy as np

### BASE FUNCTIONS

    # For string manipulations
    
def del_str(string,num):
    """Dedicated function to delete the character in position *num* in a string"""
    string=list(string)
    del(string[num])
    res=""
    for elt in string:
        res+=elt
    return res
    
def find_dash(string):
    """returns the positions of all dashes in a string"""
    res=[]
    string = list(string)
    for i,elt in enumerate(string):
        if elt =="-":
            res.append(i)
    return res
    
def insert_str(string,num,char):
    """dedicated function to insert the character *char* at position *num* in *string*"""
    string=list(string)
    string.insert(num,char)
    res=""
    for elt in string:
        res+=elt
    return res

    # For sequence dictionnary manipulations

def insert_gap_seq_dict(dict_seq,pos):
    """
    Dedicated function to add a gap in a sequence dictionnary at the position *pos*.
    """
    new_dict_seq={}
    first=True
    for in_pos, value in dict_seq.items():
        if in_pos >= pos:
            if first:
                new_dict_seq[pos]="-"
                first=False
            new_dict_seq[in_pos+1]=value
        elif in_pos < pos:
            new_dict_seq[in_pos]=value
            
    return new_dict_seq

def dict_seq_translation(dict_seq,trans_amount):
    """
    Dedicated function to translate a sequence dictionnary by *trans_amount* positions.
    """
    new_dict_seq={}
    for elt in dict_seq.items():
        new_dict_seq[trans_amount+elt[0]]=elt[1]
    
    return new_dict_seq

def dict_seq_reagglomerate(dict_seq):
    """
    Dedicated function to retransform a sequence dictionnary into an actual sequence string.
    """
    seq=""
    for elt in dict_seq.values():
        seq+=elt
    return seq

def dict_seq_cluster_finder(dict_seq,direction="Both"):
    """
    Function used to find clusters of parentheses in a sequence_dictionnary.
    """
    if direction=="Right":
        par=")"
    elif direction=="Left":
        par="("
    else:
        par="()"
        
    cluster_list_conca=[]
    for elt in dict_seq.items():
        if elt[1] in par:
            cluster_list_conca.append(elt[0])

    
    return cluster_list_conca

### PATTERN ALIGNMENT

def middle_aligning(dict_seq1,dict_seq2, diff1, diff2, mid_g1, mid_d1, mid_g2, mid_d2):
    """
    Function used to align the middle of patterns.
    
    Used in inside-out aligning.
    """
    trans_amnt=abs(mid_g1-mid_g2) 
    
    if diff1==diff2:
        #only translation operation
        if mid_g1<mid_g2:
            #middle of 1 starts before middle of 2, so 1 is to be translated
            dict_seq1=dict_seq_translation(dict_seq1, trans_amnt)
            mid_g1+=trans_amnt
            mid_d1+=trans_amnt
        else:
            #middle of 2 starts before middle of 1, so 2 is to be translated
            dict_seq2=dict_seq_translation(dict_seq2, trans_amnt)
            mid_g2+=trans_amnt
            mid_d2+=trans_amnt

    elif diff1>diff2:
        #middle of pat1 > middle of pat2
        for i in range(diff1-diff2):
            dict_seq2=insert_gap_seq_dict(dict_seq2, mid_d2)
        
        mid_d2=mid_d2+(diff1-diff2)
        
        if mid_g1<mid_g2:
            #middle of 1 starts before middle of 2, so 1 is to be translated
            dict_seq1=dict_seq_translation(dict_seq1, trans_amnt)
            mid_g1+=trans_amnt
            mid_d1+=trans_amnt
        else:
            #middle of 2 starts before middle of 1, so 2 is to be translated
            dict_seq2=dict_seq_translation(dict_seq2, trans_amnt)
            mid_g2+=trans_amnt
            mid_d2+=trans_amnt
    
    elif diff1<diff2:
        #middle of pat2 > middle of pat1
        for i in range(diff2-diff1):
            dict_seq1=insert_gap_seq_dict(dict_seq1, mid_d1)
        
        mid_d1=mid_d1+(diff2-diff1)
        
        if mid_g1<mid_g2:
            #middle of 1 starts before middle of 2, so 1 is to be translated
            dict_seq1=dict_seq_translation(dict_seq1, trans_amnt)
            mid_g1+=trans_amnt
            mid_d1+=trans_amnt
        else:
            #middle of 2 starts before middle of 1, so 2 is to be translated
            dict_seq2=dict_seq_translation(dict_seq2, trans_amnt)
            mid_g2+=trans_amnt
            mid_d2+=trans_amnt
        
    return dict_seq1, dict_seq2, mid_g1, mid_d1, mid_g2, mid_d2

def propagation_alignment(dict_tba1, dict_tba2, direction):
    """
    Aligns almost optimally two sequence_dictionnary slices from left to right or right to left.
    """
    cluster_1 = dict_seq_cluster_finder(dict_tba1, direction)
    cluster_2 = dict_seq_cluster_finder(dict_tba2, direction)

    Pdiff = len(cluster_1)-len(cluster_2)
    
    if direction=="Right" or Pdiff==0:
        #aligning the right dictionnaries.
        if Pdiff>0:
            #more pairings in 1.
            for i, elt2 in enumerate(cluster_2):
                if elt2 != cluster_1[i]:
                    local_diff = elt2-cluster_1[i]
                    if local_diff > 0:
                        #place local_diff gaps in 1 at pos cluster_1[i].
                        for k in range(local_diff):
                            dict_tba1=insert_gap_seq_dict(dict_tba1, cluster_1[i])
                            
                        #update cluster_1
                        for j in range(i, len(cluster_1)):
                            cluster_1[j]+=local_diff
                    else:
                        #place abs(local_diff) gaps in 2 at pos elt2.
                        for k in range(abs(local_diff)):
                            dict_tba2=insert_gap_seq_dict(dict_tba2, elt2)
                            
                        #update cluster_2
                        for j in range(i, len(cluster_2)):
                            cluster_2[j]+=abs(local_diff)
                        
        else:
            #more pairings in 2 or equal amount.
            for i, elt1 in enumerate(cluster_1):
                if elt1 != cluster_2[i]:
                    local_diff=cluster_2[i]-elt1
                    if local_diff > 0:
                        #place local_diff gaps in 1 at pos elt1.
                        for k in range(local_diff):
                            dict_tba1=insert_gap_seq_dict(dict_tba1, elt1)
                            
                        #update cluster_1
                        for j in range(i, len(cluster_1)):
                            cluster_1[j]+=local_diff
                    else:
                        #place abs(local_diff) gaps in 2 at pos cluster_2[i].
                        for k in range(abs(local_diff)):
                            dict_tba2=insert_gap_seq_dict(dict_tba2, cluster_2[i])
                            
                        #update cluster_2
                        for j in range(i, len(cluster_2)):
                            cluster_2[j]+=abs(local_diff)
                        
    elif direction=="Left" and Pdiff!=0:
        #aligning the left dictionnaries by starting from the right.

        if Pdiff>0:
            #more pairings in 1.
            for i in range(len(cluster_2)-1,-1,-1):
                if cluster_2[i] != cluster_1[i+Pdiff]:
                    local_diff = cluster_2[i]-cluster_1[i+Pdiff]
                    if local_diff > 0:

                        for k in range(local_diff):
                            dict_tba2=insert_gap_seq_dict(dict_tba2, cluster_2[i]+1)
                            

                        for j in range(i+1,-1,-1):
                            cluster_2[j]-=local_diff
                    else:

                        for k in range(abs(local_diff)):
                            dict_tba1=insert_gap_seq_dict(dict_tba1, cluster_1[i+Pdiff]+1)

                        for j in range(i+1,-1,-1):
                            cluster_1[j]-=abs(local_diff)
        else:
            #more pairings in 2 or equal amount.
            for i in range(len(cluster_1)-1,-1,-1):

                if cluster_2[i+abs(Pdiff)] != cluster_1[i]:
                    
                    local_diff = cluster_2[i+abs(Pdiff)]-cluster_1[i]

                    if local_diff > 0:

                        for k in range(local_diff):
                            dict_tba2=insert_gap_seq_dict(dict_tba2, cluster_2[i+abs(Pdiff)]+1)
                            
                        for j in range(i+1,-1,-1):
                            cluster_2[j]-=local_diff
                    else:

                        for k in range(abs(local_diff)):
                            dict_tba1=insert_gap_seq_dict(dict_tba1, cluster_1[i]+1)

                        for j in range(i+1,-1,-1):
                            cluster_1[j]-=abs(local_diff)


    #account for the rest of pairings by placing gaps.
    
    main_diff = len(dict_tba1)-len(dict_tba2)

    finish1=list(dict_tba1.keys())[-1]
    finish2=list(dict_tba2.keys())[-1]
    start1=list(dict_tba1.keys())[0]
    start2=list(dict_tba2.keys())[0]

    if main_diff==0:
        return dict_tba1, dict_tba2
    elif direction == "Right":
        if main_diff > 0:
            #gaps to be placed in 2
            for i in range(main_diff):
                dict_tba2[finish2+i+1]="-"
    
        elif main_diff < 0:
            #gaps to be placed in 1
            for i in range(abs(main_diff)):
                dict_tba1[finish1+i+1]="-"
    elif direction=="Left":
        if main_diff > 0:
            #gaps to be placed in 2
            for i in range(main_diff):
                dict_tba2=insert_gap_seq_dict(dict_tba2,start2)
    
        elif main_diff < 0:
            #gaps to be placed in 1
            for i in range(abs(main_diff)):
                dict_tba1=insert_gap_seq_dict(dict_tba1,start1)
    
    return dict_tba1, dict_tba2



def inside_out_pat_alignment(pat1, pat2):
    """
    Aligns matched patterns with an inside out method:
    
    First find and aligns the middle of both patterns.
    
    Then aligns the left of both patterns.
    
    In the end, aligns the right of both patterns.
    """

    #creating base dictionnary to manipulate.

    par="()"   
    dict_par1={}
    dict_seq1={}
    dict_par2={}
    dict_seq2={}
    for i,elt in enumerate(pat1.sequence):
        if elt in par:
            dict_par1[i]=elt
        dict_seq1[i]=elt
    
    for i,elt in enumerate(pat2.sequence):
        if elt in par:
            dict_par2[i]=elt
        dict_seq2[i]=elt
    
    # determining the middle of both patterns.
    
    L1 = list(dict_par1.items())
    L2 = list(dict_par2.items())
    
    i=0
    while L1[i+1][1]!=")":
        i+=1
    mid_g1 = L1[i][0]
    mid_d1 = L1[i+1][0]
    
    i=0
    while L2[i+1][1]!=")":
        i+=1
    mid_g2 = L2[i][0]
    mid_d2 = L2[i+1][0]
    
    #aligning the middle
    diff1=mid_d1-mid_g1
    diff2=mid_d2-mid_g2

    dict_seq1, dict_seq2, mid_g1, mid_d1, mid_g2, mid_d2 = middle_aligning(dict_seq1, dict_seq2, diff1, diff2, mid_g1, mid_d1, mid_g2, mid_d2)

    #slicing the patterns in 3 parts: Left, Middle, Right
    dict_tba1_L = {}
    dict_tba1_R = {}
    dict_mid1={}
    for elt in dict_seq1.items():
        if elt[0]<=mid_g1:
            dict_tba1_L[elt[0]]=elt[1]
        elif elt[0]>=mid_d1:
            dict_tba1_R[elt[0]]=elt[1]
        else:
            dict_mid1[elt[0]]=elt[1]
    
    dict_tba2_L = {}
    dict_tba2_R = {}
    dict_mid2 = {}
    for elt in dict_seq2.items():
        if elt[0]<=mid_g2:
            dict_tba2_L[elt[0]]=elt[1]
        elif elt[0]>=mid_d2:
            dict_tba2_R[elt[0]]=elt[1]
        else:
            dict_mid2[elt[0]]=elt[1]
    
    #aligning the left dictionnaries
    
    Pdiff_L = ( (dict_seq_reagglomerate(dict_tba1_L).count("(")+dict_seq_reagglomerate(dict_tba1_L).count(")")) 
               - (dict_seq_reagglomerate(dict_tba2_L).count("(")+dict_seq_reagglomerate(dict_tba2_L).count(")")) )
    
    
    if Pdiff_L == 0:
        start1=list(dict_tba1_L)[0]
        dict_tba1_L_translated=dict_seq_translation(dict_tba1_L,-start1)
        start2=list(dict_tba2_L)[0]
        dict_tba2_L_translated=dict_seq_translation(dict_tba2_L,-start2)
            
        #dict1_L, dict2_L = equal_propagation_alignment(dict_tba1_L_translated, dict_tba2_L_translated)
        dict1_L, dict2_L = propagation_alignment(dict_tba1_L_translated, dict_tba2_L_translated, "Left")
    else:
        dict1_L, dict2_L = propagation_alignment(dict_tba1_L, dict_tba2_L, "Left")

        start1=list(dict1_L.keys())[0]
        start2=list(dict2_L.keys())[0]
        
        dict1_L=dict_seq_translation(dict1_L,-start1)
        dict2_L=dict_seq_translation(dict2_L,-start2)
        
    # Now, translate the middle and right dictionnaries
    trans1 = list(dict1_L.values()).count("-") - start1
    trans2 = list(dict2_L.values()).count("-") - start2
    
    dict_tba1_R = dict_seq_translation(dict_tba1_R,trans1)
    dict_tba2_R = dict_seq_translation(dict_tba2_R,trans2)
    
    dict1_M = dict_seq_translation(dict_mid1, trans1)
    dict2_M = dict_seq_translation(dict_mid2, trans2)
    
    #aligning the right dictionnaries

    dict1_R, dict2_R = propagation_alignment(dict_tba1_R, dict_tba2_R, "Right")
    
    #reagglomerating.
    seq1_L=dict_seq_reagglomerate(dict1_L)
    seq1_M=dict_seq_reagglomerate(dict1_M)
    seq1_R=dict_seq_reagglomerate(dict1_R)
    
    seq1=seq1_L+seq1_M+seq1_R
    
    seq2_L=dict_seq_reagglomerate(dict2_L)
    seq2_M=dict_seq_reagglomerate(dict2_M)
    seq2_R=dict_seq_reagglomerate(dict2_R)
    
    seq2=seq2_L+seq2_M+seq2_R

    return seq1, seq2

def pattern_alignment(struct1, struct2, pat1, pat2, order1, order2, verbose=False):
    """
    Used to align two patterns with inside-out alignment and update the states of Structure and Pattern objects.
    """
    
    if verbose:
        print("Aligning:",pat1.pattern_nb,"with",pat2.pattern_nb)
    
    if pat1.sequence==pat2.sequence:
        if verbose:
            print("Same pattern")
        pat1.aligned(pat2,pat2.sequence)
        pat2.aligned(pat1,pat1.sequence)
        
    elif pat1.start == pat2.start and pat1.length == pat2.length:
        if verbose:
            print("Already the same length.")
        pat1.aligned(pat2,pat2.sequence)
        pat2.aligned(pat1,pat1.sequence)
    else:
        seq1,seq2=inside_out_pat_alignment(pat1, pat2)
        
        pat1.aligned(pat2,seq1)
        pat2.aligned(pat1,seq2)
        
        added_gaps1=pat1.alignedsequence.count("-")
        added_gaps2=pat2.alignedsequence.count("-")

        add_gaps(pat1.nb, added_gaps1, order1)
        struct1.length+=added_gaps1
        
        add_gaps(pat2.nb, added_gaps2, order2)
        struct2.length+=added_gaps2

### STRUCTURE INITIALIZING

def subdiv_finder(sequence,subdiv_param):
    """
    Finds and marks with hashtags the pairings that compose overdivisions.
    """
    par="()"
    #void=".-"
    dict_par={}
    dict_seq={}
    for i,elt in enumerate(sequence):
        if elt in par:
            dict_par[i]=elt
        dict_seq[i]=elt
        
    new_dict={}
    subdiv={}
    while dict_par != {}:
        index_list=list(dict_par.keys())
        i=0
        while dict_par[index_list[i]]==dict_par[index_list[i+1]]:
            i+=1
        if not((index_list[i+1] - index_list[i]) > int(len(sequence)/subdiv_param)):
            new_dict[index_list[i]]=dict_par[index_list[i]]
            new_dict[index_list[i+1]]=dict_par[index_list[i+1]]
        else:
            subdiv[index_list[i]]=dict_par[index_list[i]]
            subdiv[index_list[i+1]]=dict_par[index_list[i+1]]
            
        del(dict_par[index_list[i]])
        del(dict_par[index_list[i+1]])
        
    for elt in subdiv.keys():
        dict_seq[elt]="#"
    new_seq=''
    for elt in dict_seq.values():
        new_seq+=elt
    
    
    order=[]
    for elt in subdiv.items():
        order.append(elt)
    
    def get_index(pair):
        return pair[0]
    
    order=sorted(order, key=lambda pair : get_index(pair))
    
    subdiv_fin={}
    for i,pairs in enumerate(order):
        subdiv_fin[i]=pairs[1]
        
    return subdiv_fin,new_seq

def slicer(sequence):
    """
    Slices Structure objects into Patterns and Separators.
    
    Used in the Structure class constructor.
    """
    pat=[]
    sep=[]
    par="()"
    dict_par={}
    dict_seq={}
    for i,elt in enumerate(sequence):
        if elt in par:
            dict_par[i]=elt
        dict_seq[i]=elt
    open_db = 0
    start = min(dict_par.keys())
    finish=0
    new=True
    ranges=[]
    for elt in dict_par.keys():
        if dict_par[elt]=="(":
            if new:
                start=elt
                new=False
            open_db+=1
        elif dict_par[elt]==")":
            open_db-=1
        if open_db==0:
            finish=elt
            ranges.append((start,finish+1))
            new=True
    
    full_ranges=[]
    for tup in ranges:
        full_ranges.append(("PAT",[i for i in range(tup[0],tup[1])]))
    
    if sequence[0]!="(":
        separator_start = True
    else:
        separator_start = False
        
    if sequence[-1]!=")":
        separator_end = True
    else:
        separator_end = False
        
    dp_full_ranges=full_ranges[:]
    insert=0
    if separator_start:
        dp_full_ranges.insert(0,("SEP",[j for j in range(full_ranges[0][1][0])]))
        insert+=1
    else:
        dp_full_ranges.insert(0,("SEP",[]))
        insert+=1
    for i in range(1,len(full_ranges)):
        dp_full_ranges.insert(i+insert,("SEP",[j for j in range(full_ranges[i-1][1][-1]+1,full_ranges[i][1][0])]))
        insert+=1
    if separator_end:
        dp_full_ranges.append(("SEP",[j for j in range(full_ranges[-1][1][-1]+1,len(sequence))]))
    else:
        dp_full_ranges.append(("SEP",[]))
    
    order=0
    pat_num=0
    for i,elt in enumerate(dp_full_ranges):
        if elt[0]=="SEP":
            new_seq=''
            if not(elt[1] == []):
                for nbs in elt[1]:
                    new_seq+=dict_seq[nbs]
                sep.append(Separator(new_seq,(elt[1][0],elt[1][-1]),order))
            else:
                if order==0:
                    sep.append(Separator(new_seq,(-1,0),order))
                else:
                    sep.append(Separator(new_seq,(-1,dp_full_ranges[i-1][1][-1]+1),order))
            order+=1
        if elt[0]=="PAT":
            new_seq=''
            for nbs in elt[1]:
                new_seq+=dict_seq[nbs]
            pat.append(Pattern(new_seq,(elt[1][0],elt[1][-1]),order,pat_num))
            pat_num+=1
            order+=1
    
    return sep,pat

### CLASSES

class Structure():
    """
    Class representing an aligned or not dotbracket sequence.
    
    It is formed of two lists : one with its patterns and one with its separators.
    """
    def __init__(self, sequence, ident=None, fam=None, AGU=None):
        self.raw=sequence
        self.sequence=AGU
        self.subdiv,self.raw_nosubdiv=subdiv_finder(sequence, 2)
        sep,pat=slicer(self.raw_nosubdiv)
        self.separators=sep
        self.patterns=pat
        self.length=len(sequence)
        self.pattern_nb=len(pat)
        self.separator_nb=len(sep)
        self.isaligned=False
        self.id=ident
        self.family=fam
        self.alignedsequence=""
        self.alignedwith=None
    
    def reset(self, sequence):
        self.raw=sequence
        self.subdiv,self.raw_nosubdiv=subdiv_finder(sequence, 2)
        sep,pat=slicer(self.raw_nosubdiv)
        self.separators=sep
        self.patterns=pat
        self.length=len(sequence)
        self.pattern_nb=len(pat)
        self.separator_nb=len(sep)
        self.isaligned=False
        self.alignedsequence=""
        self.alignedwith=None
    
    def __str__(self):
        tbp="Type: Structure\n"
        tbp+="ID: "+str(self.id)+"\n"
        tbp+="Length: "+str(self.length)+"\n"
        tbp+="Raw Sequence: "+self.raw+"\n"
        tbp+="Subdiv sequence: "+str(self.raw_nosubdiv)+"\n"
        tbp+="Subdiv: "+str(self.subdiv)+"\n"
        tbp+="\nSeparators: \n"
        for i,elt in enumerate(self.separators):
            tbp+="\nSeparator number "+str(i+1)+":\n"+str(elt)+"\n"
        tbp+="\nPatterns: \n"
        for i,elt in enumerate(self.patterns):
            tbp+="\nPattern number "+str(i+1)+":\n"+str(elt)+"\n"

        if self.isaligned:
            tbp+="Aligned sequence: "+self.alignedsequence+"\n"
            tbp+="Aligned with:"+self.alignedwith.alignedsequence+"\n"
        else:
            tbp+="Not yet aligned\n"
        return tbp

    def aligned(self):
        self.isaligned=True
        for pattern in self.patterns:
            if not pattern.isaligned:
                self.isaligned=False
        if not self.isaligned:
            print("Structure not yet aligned")
        

    def __eq__(self, other):
        if isinstance(other, Structure):
            return other.raw == self.raw
        
    def order_list(self):
        
        def get_order(motif):
            return motif.nb
        
        non_ordered=[]
        for pattern in self.patterns:
                non_ordered.append(pattern)
        for separator in self.separators:
            non_ordered.append(separator)
            
        ordered = sorted(non_ordered, key=lambda pt : get_order(pt))
        
        return ordered
    
    def reagglomerate(self):
        if self.isaligned:
            
            for pattern in self.patterns:
                if not pattern.isaligned:
                    print("The patterns are not yet aligned, returning raw sequence")
                    return self.raw
            
            ordered=self.order_list()
            
            seqint=""
            for patsep in ordered:
                if isinstance(patsep,Separator):
                    seqint+=patsep.sequence
                else:
                    seqint+=patsep.alignedsequence
            
            seq_dic={}
            cpt=0
            for i,elt in enumerate(seqint):
                if elt == "#":
                    seq_dic[i]=self.subdiv[cpt]
                    cpt+=1
                else:
                    seq_dic[i]=elt
            
            seq=''
            for elt in seq_dic.values():
                seq+=elt
            
            return seq
        
        else:
            print("The structure is not yet aligned, returning raw sequence")
            return self.raw

class Separator():
    """
    Class used to define a structure.
    A separator is represented by an array of points separating two patterns.
    
    Example: 
        In the structure `(((...)))...(.(((...))))` the three dots in the middle forms the separator.
    """
    
    def __init__(self,raw,raw_range,order):
        self.nb=order
        self.start=raw_range[0]
        self.finish=raw_range[1]
        if raw_range[0] == -1 :
            self.length = 0
        else:
            self.length=raw_range[1]-raw_range[0]+1
            
        self.sequence=raw
        
        self.subdiv_index=raw.count('#')
        
        if self.subdiv_index!=0:
            self.subdiv_accounted=False
            
    def __eq__(self,other):
        if isinstance(other,Separator):
            return self.length==other.length
    
    def __str__(self):
        tbp="Type: Separator\n"
        tbp+="Length: "+str(self.length)+"\n"
        tbp+="Order: "+str(self.nb)+"\n"
        tbp+="Sequence: "+self.sequence+"\n"
        tbp+="Start: "+str(self.start)+"\n"
        tbp+="Finish: "+str(self.finish)+"\n"
        tbp+="Subdiv index: "+str(self.subdiv_index)+"\n"
        return tbp
    
    def compare(self,other):
        if isinstance(other,Separator):
            return abs(self.length-other.length)

class Pattern():
    """
    Class used to define a structure.
    A pattern is represented by a smaller structure with a single opening and closing sequence.
    
    Example: `(((((..))..)))` is a pattern but `(..)(((.)))` isn't.
    """
    def __init__(self,raw,raw_range,order,pat_num):
        self.nb=order
        self.pattern_nb=pat_num
        self.start=raw_range[0]
        self.finish=raw_range[1]
        self.length=raw_range[1]-raw_range[0]+1
        self.sequence=raw
        
        self.subdiv_index=raw.count("#")
        
        self.isaligned=False
        self.alignedwith=None
        self.alignedsequence=""

    def __eq__(self,other):
        if isinstance(other,Pattern):
            return self.sequence==other.sequence
    
    def am_distance(self,other):
        if isinstance(other,Pattern):
            return AF.compute_distance_clustering(AF.SecondaryStructure(self.sequence),AF.SecondaryStructure(other.sequence), "cityblock", "slow")

    def aligned(self,acc_pat,new_seq):
        self.isaligned=True
        if acc_pat==None:
            self.alignedwith=EmptyPattern
        else:
            self.alignedwith=acc_pat
        self.alignedsequence=new_seq
        self.length=len(new_seq)
        self.finish+=new_seq.count('-')
        
    def __str__(self):
        tbp="Type: Pattern\n"
        tbp+="Length: "+str(self.length)+"\n"
        tbp+="Order: "+str(self.nb)+"\n"
        tbp+="Start: "+str(self.start)+"\n"
        tbp+="Finish: "+str(self.finish)+"\n"
        tbp+="Starting Sequence: "+self.sequence+"\n"
        tbp+="Subdiv index: "+str(self.subdiv_index)+"\n"
        if self.isaligned:
            tbp+="Aligned sequence: "+self.alignedsequence+"\n"
            tbp+="Aligned with: "+self.alignedwith.alignedsequence+"\n"
        else:
            tbp+="Not yet aligned\n"
        return tbp

EmptyPattern=Pattern('',[-1,-1],-1,-1)

### BASE SEPARATOR GAPS FUNCTIONS=

def add_gaps(current_order,added_gaps,order_list):
    """
    Updates the start and finish positions of further separators and patterns.
    """
    for elt in order_list:
        if elt.nb > current_order:
            #Test for an update of an empty separator.
            if elt.start==-1:
                elt.finish+=added_gaps
            else:
                elt.start+=added_gaps
                elt.finish+=added_gaps

def sep_gap_adder(struct,nb_gaps,sep,ordered):
    """
    Adds `nb_gaps` gaps at the end of `sep` and updates the other pat and sep with `ordered` and `struct`
    """
    
    #updating positions for an empty separator or a normal separator.
    if sep.start==-1:
        sep.start=sep.finish
        sep.finish+=nb_gaps
    else:
        sep.finish+=nb_gaps
    sep.sequence+=nb_gaps*'-'
    add_gaps(sep.nb, nb_gaps, ordered)
    struct.length+=nb_gaps

### SEPARATOR ALIGNMENT FUNCTIONs

def overdivision_compensating(struct1, struct2, ordered1, ordered2, matching):
    """
    Accounting for overdivision and aligning separators where hashtags are present.
    
    We are only looking at separators present right before or after matched patterns.
    
    If there are hashtags in both separators compared, we align them with eachother.
    """
    
    #creating necessary separator matching 
    sep_matching=[]
    
    for pairs in matching:
        if ordered1[pairs[0].nb-1].subdiv_index != 0 and ordered2[pairs[1].nb-1].subdiv_index !=0:
            if not(ordered1[pairs[0].nb-1].subdiv_accounted) and not(ordered2[pairs[1].nb-1].subdiv_accounted):
                ordered1[pairs[0].nb-1].subdiv_accounted = True
                ordered2[pairs[1].nb-1].subdiv_accounted = True
                sep_matching.append([ordered1[pairs[0].nb-1],ordered2[pairs[1].nb-1]])
            
        if ordered1[pairs[0].nb+1].subdiv_index != 0 and ordered2[pairs[1].nb+1].subdiv_index !=0:
            if not(ordered1[pairs[0].nb+1].subdiv_accounted) and not(ordered2[pairs[1].nb+1].subdiv_accounted):
                ordered1[pairs[0].nb+1].subdiv_accounted = True
                ordered2[pairs[1].nb+1].subdiv_accounted = True
                sep_matching.append([ordered1[pairs[0].nb+1],ordered2[pairs[1].nb+1]])

    for sep_pairs in sep_matching:
        #making sequence dictionnaries
        dict_seq1={}
        for i,elt in enumerate(sep_pairs[0].sequence):
            dict_seq1[i]=elt
        
        dict_seq2={}
        for i,elt in enumerate(sep_pairs[1].sequence):
            dict_seq2[i]=elt
        
        #determining the start of both subdiv
        start1=0
        while dict_seq1[start1]!="#":
            start1+=1
            
        start2=0
        while dict_seq2[start2]!="#":
            start2+=1
        
        
        #aligning if subdiv starts at different places in the separator
        if start1!=start2:
            
            diff=start1-start2
            
            if diff > 0:
                #start1>start2, # starts after in sep1; placing gaps in sep2.
                for i in range(abs(diff)):
                    dict_seq2 = insert_gap_seq_dict(dict_seq2,start2)
            elif diff <0:
                #start2>start1, # starts after in sep2; placing gaps in sep1.
                for i in range(abs(diff)):
                    dict_seq1 = insert_gap_seq_dict(dict_seq1,start1)
            
            sep_pairs[0].sequence=dict_seq_reagglomerate(dict_seq1)
            sep_pairs[1].sequence=dict_seq_reagglomerate(dict_seq2)
            
            nb_gaps1=sep_pairs[0].sequence.count('-')
            nb_gaps2=sep_pairs[1].sequence.count('-')
            
            sep_pairs[0].length+=nb_gaps1
            sep_pairs[1].length+=nb_gaps2
            
            sep_pairs[0].finish+=nb_gaps1
            sep_pairs[1].finish+=nb_gaps2

            add_gaps(sep_pairs[0].nb,nb_gaps1,ordered1)
            add_gaps(sep_pairs[1].nb,nb_gaps2,ordered2)

            struct1.length+=nb_gaps1
            struct2.length+=nb_gaps2

def sep_gap_inserter(struct1, struct2, matching, ordered1, ordered2, main_diff):
    """
    Function used to insert gaps in separators where it is necessary in regards to:
        - the biggest structure
        - the pattern matching
    """
    
    for elt in matching:
        
        if struct1.length-struct2.length > 0:
            bigger = 1
        else:
            bigger = 2
        #struct1 ~ elt[0] - struct2 ~ elt[1]
        if not elt[0].start == elt[1].start:
            diff=elt[0].start-elt[1].start
            if diff < 0:
                #start pat1 < start pat2
                #pat1 starts before pat2
                #input diff gaps before pat1 and after pat2 if main_diff=0 or bigger=1 else only before pat1.

                if main_diff==0:
                    sep_gap_adder(struct1, abs(diff), ordered1[elt[0].nb-1],ordered1)
                    sep_gap_adder(struct2, abs(diff), ordered2[elt[1].nb+1],ordered2)
                elif bigger == 2:
                    sep_gap_adder(struct1, abs(diff), ordered1[elt[0].nb-1],ordered1)
                elif bigger == 1:
                    sep_gap_adder(struct1, abs(diff), ordered1[elt[0].nb-1],ordered1)
                    sep_gap_adder(struct2, abs(diff), ordered2[elt[1].nb+1],ordered2)
            else:
                #start pat1 > start pat2
                #pat2 starts before pat1
                #input diff gaps after pat1 and before pat2 if main_diff=0 or bigger=2 else only before pat2

                if main_diff==0:
                    sep_gap_adder(struct2, abs(diff), ordered2[elt[1].nb-1],ordered2)
                    sep_gap_adder(struct1, abs(diff), ordered1[elt[0].nb+1],ordered1)
                elif bigger==2:
                    sep_gap_adder(struct2, abs(diff), ordered2[elt[1].nb-1],ordered2)
                    sep_gap_adder(struct1, abs(diff), ordered1[elt[0].nb+1],ordered1)
                elif bigger==1:
                    sep_gap_adder(struct2, abs(diff), ordered2[elt[1].nb-1],ordered2)
                    
        main_diff=abs(struct1.length-struct2.length)
              
    #Adding the last gaps at the end of the structures if main_diff is still positive.

    if struct1.length-struct2.length > 0:
        bigger = 1
    else:
        bigger = 2
    
    if main_diff>0:
        if bigger == 1:
            last_sep=ordered2[-1]
            added_gaps2=main_diff
            last_sep.sequence+=added_gaps2*'-'
            add_gaps(last_sep.nb, added_gaps2, ordered2)
            struct2.length+=added_gaps2
        elif bigger==2:
            last_sep=ordered1[-1]
            added_gaps1=main_diff
            last_sep.sequence+=added_gaps1*'-'
            add_gaps(last_sep.nb, added_gaps1, ordered1)
            struct1.length+=added_gaps1

def separator_compensating(struct1, struct2, matching):
    """
    Used after pattern aligning to add gaps in separators to match the final length and returning the aligned structures.
    
    """
    ordered1=struct1.order_list()
    ordered2=struct2.order_list()
    
    if struct1.length==struct2.length:
        #the structures have the same size
        sep_gap_inserter(struct1, struct2, matching, ordered1, ordered2,0)
    
    else :
        main_diff=abs(struct2.length-struct1.length)
        sep_gap_inserter(struct1, struct2, matching, ordered1, ordered2, main_diff)

### MATCHING FUNCTION

def pair_pat_score_pass1(pat1,pat2):
    """
    Gives a pairing score for first pass pattern recognition based on the aptamat distance, the number of pairings and the length difference.
    """
    apta_dist=AF.compute_distance_clustering(AF.SecondaryStructure(pat1.sequence),AF.SecondaryStructure(pat2.sequence), "cityblock", "slow")
    
    pat1_sc = (pat1.sequence.count("(")+pat1.sequence.count(")"))
    pat2_sc = (pat2.sequence.count("(")+pat2.sequence.count(")"))

    return abs(pat1_sc - pat2_sc) + apta_dist + abs(pat1.length - pat2.length)

def pair_pat_score_pass2(pat1,pat2):
    """
    Gives a pairing score for second pass pattern recognition based on the aptamat distance and the number of pairings.
    """
    apta_dist=AF.compute_distance_clustering(AF.SecondaryStructure(pat1.sequence),AF.SecondaryStructure(pat2.sequence), "cityblock", "slow")
    
    pat1_sc = (pat1.sequence.count("(")+pat1.sequence.count(")"))
    pat2_sc = (pat2.sequence.count("(")+pat2.sequence.count(")"))

    return abs(pat1_sc - pat2_sc) + apta_dist

def matching_test(matching):
    pair_order=[matching[0][0].nb,matching[0][1].nb]
    
    for i in range(1,len(matching)):
        if matching[i][0].nb<pair_order[0] or matching[i][1].nb < pair_order[1]:
            return False
        else:
            pair_order[0]=matching[i][0].nb
            pair_order[1]=matching[i][1].nb
            
    return True

def naive_matching(order1, order2):
    patterns1=[]
    patterns2=[]
    for elt in order1:
        if isinstance(elt, Pattern):
            patterns1.append(elt)
            
    for elt in order2:
        if isinstance(elt, Pattern):
            patterns2.append(elt)
    
    matching=[]
    if len(patterns1)<=len(patterns2):
        for i,elt in enumerate(patterns1):
            matching.append([elt,patterns2[i]])
    elif len(patterns2) < len(patterns1):
        for i, elt in enumerate(patterns2):
            matching.append([patterns1[i],elt])
    
    return matching

def matching_finder(struct1, struct2, verbose=False):
    """
    Used to pair up patterns.
    
    MATCHING METHOD:
        - Finds the best pairings based on a pairing score calculated in another function.
        
    - Does three passes
        - 1 - with score with length difference.
        - 2 - with score without length difference.
        - 3 - In order.
        
    - Keeps the best matching.
    """
    order1=struct1.order_list()
    order2=struct2.order_list()
    #calculating best match based on pair by pair score.
    matching_pass1=[]
    matching_pass2=[]
    matching_pass3=[]
    
    pairing_matrix_pass1=np.empty((len(struct1.patterns),len(struct2.patterns)),dtype=tuple)
    pairing_matrix_score_pass1=np.empty((len(struct1.patterns),len(struct2.patterns)),dtype=float)
    
    pairing_matrix_pass2=np.empty((len(struct1.patterns),len(struct2.patterns)),dtype=tuple)
    pairing_matrix_score_pass2=np.empty((len(struct1.patterns),len(struct2.patterns)),dtype=float)
    # determining the matrix of paired scores:
    if verbose:
        print("Calculating the three passes.")
    for i,pat1 in enumerate(struct1.patterns):
        for j,pat2 in enumerate(struct2.patterns):
            pairing_matrix_pass1[i,j]=(pair_pat_score_pass1(pat1, pat2),pat1,pat2)
            pairing_matrix_score_pass1[i,j]=pair_pat_score_pass1(pat1, pat2)
            
            pairing_matrix_pass2[i,j]=(pair_pat_score_pass2(pat1, pat2),pat1,pat2)
            pairing_matrix_score_pass2[i,j]=pair_pat_score_pass2(pat1, pat2)
            
    max_pairings = min(len(struct1.patterns),len(struct2.patterns))
    
    for i in range(max_pairings):
        #PASS1
        argkeep1=np.unravel_index(np.argmin(pairing_matrix_score_pass1),pairing_matrix_score_pass1.shape)
        keep1=pairing_matrix_pass1[argkeep1[0],argkeep1[1]]
        matching_pass1.append([keep1[1],keep1[2]])
        
        pairing_matrix_pass1=np.delete(pairing_matrix_pass1,obj=argkeep1[0],axis=0)
        pairing_matrix_pass1=np.delete(pairing_matrix_pass1,obj=argkeep1[1],axis=1)
        
        pairing_matrix_score_pass1=np.delete(pairing_matrix_score_pass1,obj=argkeep1[0],axis=0)
        pairing_matrix_score_pass1=np.delete(pairing_matrix_score_pass1,obj=argkeep1[1],axis=1)
        
        #PASS2
        argkeep2=np.unravel_index(np.argmin(pairing_matrix_score_pass2),pairing_matrix_score_pass2.shape)
        keep2=pairing_matrix_pass2[argkeep2[0],argkeep2[1]]
        matching_pass2.append([keep2[1],keep2[2]])
        
        pairing_matrix_pass2=np.delete(pairing_matrix_pass2,obj=argkeep2[0],axis=0)
        pairing_matrix_pass2=np.delete(pairing_matrix_pass2,obj=argkeep2[1],axis=1)
        
        pairing_matrix_score_pass2=np.delete(pairing_matrix_score_pass2,obj=argkeep2[0],axis=0)
        pairing_matrix_score_pass2=np.delete(pairing_matrix_score_pass2,obj=argkeep2[1],axis=1)
        
        
    #reordering matchings:
    
    def get_id_pat1(pat):
        return pat[0].nb
    
    matching_pass1=sorted(matching_pass1, key=lambda pt : get_id_pat1(pt))
    matching_pass2=sorted(matching_pass2, key=lambda pt : get_id_pat1(pt))
    
    #Determining if pass1 is a valid matching, if not, use pass2 or pass3 if pass2 is not valid.
    if not matching_test(matching_pass1):
        if verbose:
            print("Pass1 not valid")
        if not matching_test(matching_pass2):
            if verbose:
                print("Pass2 not valid, using naive ordering.")
            
            matching_pass3=naive_matching(order1, order2)
            
            full_paired=[]
            for elt in matching_pass3:
                full_paired.append(elt[0])
                full_paired.append(elt[1])
            
            for pat1 in struct1.patterns:
                if pat1 not in full_paired:
                    pat1.aligned(None,pat1.sequence)
            
            for pat2 in struct2.patterns:
                if pat2 not in full_paired:
                    pat2.aligned(None,pat2.sequence)
                    
            return matching_pass3
        else:
            if verbose:
                print("Pass2 valid, returning matching.")
            full_paired=[]
            for elt in matching_pass2:
                full_paired.append(elt[0])
                full_paired.append(elt[1])
            
            for pat1 in struct1.patterns:
                if pat1 not in full_paired:
                    pat1.aligned(None,pat1.sequence)
            
            for pat2 in struct2.patterns:
                if pat2 not in full_paired:
                    pat2.aligned(None,pat2.sequence)
                    
            return matching_pass2
    else:
        if verbose:
            print("Pass1 valid, returning matching.")
        full_paired=[]
        for elt in matching_pass1:
            full_paired.append(elt[0])
            full_paired.append(elt[1])
        
        for pat1 in struct1.patterns:
            if pat1 not in full_paired:
                pat1.aligned(None,pat1.sequence)
        
        for pat2 in struct2.patterns:
            if pat2 not in full_paired:
                pat2.aligned(None,pat2.sequence)
                
        return matching_pass1

### MAIN ALIGNMENT FUNCTION

def full_alignment(struct1, struct2, verbose=False):
    """
    
    Main function to calculate a pattern based alignment.
    
    struct1: Structure
    struct2: Structure
    
    ____________________
    METHOD:
        
        MATCHING
        
        - Testing for compatibility using a pairing score.
        
        - If the best match is not compatible with an alignment, returns the matching instead.
        
        ALIGNING
        
        - Using inside-out alignment to align each pattern with its match.
            - mark them as aligned and input each other in the `self.alignedwith` variable.
        
        - Test for overdivisions alignment possibilities and aligns the separators compatible.
        
        - When done, match the length of the two structures with gaps placed in the separators where it minimizes the distance.
        
        - Mark the structures as aligned and input each other in the `self.alignedwith` variable.
        
        - return the structures and the new distance.
    ____________________
    
    """
    if verbose:
        a=time.time()
        initial_dist=AF.compute_distance_clustering(AF.SecondaryStructure(struct1.raw),AF.SecondaryStructure(struct2.raw), "cityblock", "slow")
    
    matching=[]

    if struct1.raw == struct2.raw:
        for i,pat in enumerate(struct1.patterns):
            pat.aligned(struct2.patterns[i],pat.sequence)
            struct2.patterns[i].aligned(pat,struct2.patterns[i].sequence)
        
        struct1.alignedwith=struct2
        struct2.alignedwith=struct1
        
        struct1.aligned()
        struct2.aligned()
        
        struct1.alignedsequence = struct1.reagglomerate()
        struct2.alignedsequence = struct2.reagglomerate()
        if verbose:
            print("Same structures, returning the same sequence")
        return struct1, struct2
        
    else:
        if verbose:
            print("\nDetermining the optimal pattern match if possible \n")
        
        matching = matching_finder(struct1, struct2, verbose)
    
    if verbose:
        print("\nAligning patterns.")
    

    for elt in matching:
        order1=struct1.order_list()
        order2=struct2.order_list()
        pattern_alignment(struct1,struct2,elt[0],elt[1],order1,order2,verbose)

    if verbose:
        print("\nAccounting for overdivison")

        
    order1=struct1.order_list()
    order2=struct2.order_list()
    overdivision_compensating(struct1, struct2, order1, order2, matching)

    if verbose:
        print("\nAdding gaps in separators for length and pattern matching")

    separator_compensating(struct1, struct2, matching)

    if verbose:
        print("\nUpdating last parameters and finishing\n")
    struct1.alignedwith=struct2
    struct2.alignedwith=struct1
    
    struct1.aligned()
    struct2.aligned()
  
    struct1.alignedsequence = struct1.reagglomerate()
    struct2.alignedsequence = struct2.reagglomerate()
    
    if verbose:
        new_dist = AF.compute_distance_clustering(AF.SecondaryStructure(struct1.alignedsequence),AF.SecondaryStructure(struct2.alignedsequence), "cityblock", "slow")
        print(struct1.alignedsequence)
        print(struct2.alignedsequence)
        if initial_dist!=0:
            improvement = round((initial_dist-new_dist)/initial_dist *100,2)
        else:
            improvement = 1
        print("Improvement:",initial_dist,"->",new_dist,"| in %:",str(improvement)+"%")
        b=time.time()
        print("Time spent:",str(round(b-a,3))+"s")

#_____________________________MAIN FUNCTIONS_____________________________

def row_builder(i,struct1, structure_list):
    row=[]
    for j,struct2 in enumerate(structure_list):
        row.append((AF.compute_distance_clustering(AF.SecondaryStructure(struct1.raw),AF.SecondaryStructure(struct2.raw), "cityblock", "slow"),struct1,struct2))
    return row

def build_aptamat_matrix(structure_list,verbose):
    if verbose:
        print("Building matrix (this might take a while)")

    CORE=mp.cpu_count()
    
    pool = mp.Pool(CORE)
    matrix=[]
    matrix.append(pool.starmap(row_builder, [(i,struct1, structure_list) for i,struct1 in enumerate(structure_list)]))
    pool.terminate()
    
    return np.array(matrix[0])


def argmin(row,row_number):
    mini=row[0][0]
    argmini=0
    for i,elt in enumerate(row):
        if i!=row_number:
            if elt[0]<=mini:
                mini=elt[0]
                argmini=i
    return argmini, mini

def ensemble_aligning(matrix):
    aligned_structure_list=[]
    for i,rows in enumerate(matrix):
        argkeep,minval=argmin(rows,i)
        full_alignment(matrix[i][argkeep][1], matrix[i][argkeep][2])
        #matrix[i][argkeep][2].reset(matrix[i][argkeep][2].alignedsequence)
        aligned_structure_list.append(matrix[i][argkeep][1])

    return aligned_structure_list

def initialize_dataset(structure_file):
    """
    Function used to read from the structure file.
    """
    structure_list = []
    family = {}
    with open(structure_file, 'r') as file:
        for line in file:
            content = line.strip().split()
            
            if content:
                
                if line.startswith('FAMILY'):
                    pass
                else:
                    try:
                        family[content[0]] += 1
                    except KeyError:
                        family[content[0]] = 1
    
                if AF.Dotbracket.is_dotbracket(content[3]):
                    
                    structure_list.append((content[3],content[1].split('.')[0],content[0],content[2]))
                    
    return structure_list,family


def main():
    start=time.time()
    parser = argparse.ArgumentParser(description="AptAlign is an alignment algorithm designed around pattern recognition.\n"
                                                 "Use -fp for a file input of an ensemble of structures.\n"
                                                 "Use -s to input only two structures directly in the command line.\n"
                                                 "Use -v to toggle verbose mode.\n"
                                                 "Use -l to toggle the generation of a complete log file.\n"
                                                 "Use -pr to only generate a pattern recognition.\n")

    parser.add_argument('-v',
                        '--verbose',
                        help="Increase output verbosity.",
                        default=False,
                        action="store_true")

    parser.add_argument('-s',
                        '--structures',
                        nargs='+',
                        type=str,
                        help='two 2D structures in dotbracket notation.')

    parser.add_argument('-fp',
                        '--filepath',
                        action="store",
                        nargs='+',
                        help='Input filepath containing structures to be aligned in structures format.')
    
    parser.add_argument('-pr',
                        '--pattern_recognition',
                        help="Only generates pattern recognition",
                        action="store_true")
    
    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if args.filepath is not None:
        structure_file=""
        for elt in args.filepath:
            structure_file+=str(elt)
        
        if args.verbose:
            print("Reading dataset")
        
        structure_list, families = initialize_dataset(structure_file)
        
        if args.verbose:
            print("Slicing Structures")
        
        struct_obj_list=[]
        for struct in structure_list:
            struct_obj_list.append(Structure(struct[0],ident=str(struct[1]),fam=struct[2], AGU=struct[3]))
        
        fam_dict={}
        for family in families:
            fam_dict[family]=[]
            
        for struct in struct_obj_list:
            fam_dict[struct.family].append(struct)
        
        matrix=build_aptamat_matrix(struct_obj_list,args.verbose)
        
        aligned_structure_list=ensemble_aligning(matrix)
        
        outputfilepath=structure_file.replace(".dat","")+"_aptaligned.dat"
        
        tbw="FAMILY    PDB_chain    SEQUENCE    DOTBRACKET\n"
        
        for struct in aligned_structure_list:
            if args.verbose:
                print(struct.id,struct.alignedsequence)
            tbw+=struct.family+"    "+struct.id+"    "+struct.sequence+"    "+struct.alignedsequence+"\n"
        
        f_created=open(outputfilepath,'a')
        f_created.write(tbw)
        f_created.close()
        
        finish=time.time()
        tot=round(finish-start,2)
        print("Execution time: ",tot,"s")
        
        
        sys.exit(0)
        
        
        
    if args.structures is not None:
        struct1=Structure(args.structures[0])
        struct2=Structure(args.structures[1])
        
        full_alignment(struct1, struct2, args.verbose)
        if not args.verbose:
            print("Structure 1:")
            print(struct1.alignedsequence)
            print("Structure 2:")
            print(struct2.alignedsequence)
        else:
            print("\nDetailed results:\n")
            print("______________________________________\n")
            print("\n","First Structure","\n")
            print(struct1)
            print("______________________________________\n")
            print("\n","Second Structure","\n")
            print(struct2)
        sys.exit(0)

if __name__ == '__main__':
    main()
