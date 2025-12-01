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
import copy
from tabulate import tabulate

### DEBUGGING FUNCTIONS

def vis_align(struct1,struct2,step="NULL"):
    """
    Function used to visualize the alignment at any step in the full_alignment function.
    """
    seq1=""
    for elt in struct1.order_list():
        if isinstance(elt, Pattern):
            if elt.alignedsequence=="":
                seq1+=elt.sequence
            else:
                seq1+=elt.alignedsequence
        else:
            seq1+=elt.sequence
    seq2=""
    for elt in struct2.order_list():
        if isinstance(elt, Pattern):
            if elt.alignedsequence=="":
                seq2+=elt.sequence
            else:
                seq2+=elt.alignedsequence
        else:
            seq2+=elt.sequence
            
    print("\nSTEP:",step," | Calculated lengths: - Struct1: ",struct1.length,"- Struct2: ",struct2.length)
    print(seq1,"    Real Length of Struct1: ",len(seq1))
    print(seq2,"    Real Length of Struct2: ",len(seq2))
    if len(seq1) != struct1.length or len(seq2) != struct2.length:
        print("ERROR IN LENGTH CALCULATION AT STEP: ",step)

### ERROR HANDLING CLASSES

class MatchingError(Exception):
    pass

class PatternAlignmentError(Exception):
    pass

class PKCompensatingError(Exception):
    pass

class ODCompensatingError(Exception):
    pass

class SepCompensatingError(Exception):
    pass

class GeneralAlignmentError(Exception):
    pass

class StructBuildError(Exception):
    pass

class MultiprocessingError(Exception):
    pass

class OPKPrioError(Exception):
    pass

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
    """
    dedicated function to insert the character *char* at position *num* in *string*
    """
    string=list(string)
    string.insert(num,char)
    res=""
    for elt in string:
        res+=elt
    return res

def count_pseudo(seq):
    """
    Counts the number of pseudoknots characters in a sequence.
    """
    pseudo_char="[{<>}]"
    count=0
    for elt in seq:
        if elt in pseudo_char:
            count+=1
    return count

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

def inverser_dico(dict_seq):
    dico_tbd={}
    maxi=max(dict_seq.keys())
    mini=min(dict_seq.keys())
    for i in range(maxi,mini-1,-1):
        dico_tbd[maxi-i]=dict_seq[i]
    
    return dico_tbd

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
    Function used to find all parentheses and isolate them in a sequence_dictionnary.
    """
    if direction=="Right":
        par=")"
    elif direction=="Left":
        dict_seq=inverser_dico(dict_seq)
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
    
    #counting pseudoknots and determining start. Aligning them if possible.
    pseudo="[{<>}]"
    seq1=""
    for key,value in dict_seq1.items():
        if mid_g1 < key < mid_d1:
            seq1+=value

    seq2=""
    for key,value in dict_seq2.items():
        if mid_g2 < key < mid_d2:
            seq2+=value
    
    nb_pk1=seq1.count("[")+seq1.count("{")+seq1.count("<")+seq1.count(">")+seq1.count("}")+seq1.count("]")
    nb_pk2=seq2.count("[")+seq2.count("{")+seq2.count("<")+seq2.count(">")+seq2.count("}")+seq2.count("]")
    if nb_pk1 !=0 and nb_pk2 !=0:
        start1=0
        while seq1[start1] not in pseudo:
            start1+=1
        
        start2=0
        while seq2[start2] not in pseudo:
            start2+=1
    
        if start1<start2:
            #displace mid1 by start1 gaps at start1 to make pseudoknots match.
            for i in range(start2-start1):
                dict_seq1=insert_gap_seq_dict(dict_seq1, start1+mid_g1+1)
                  
            mid_d1=mid_d1+(start2-start1)
        
        elif start2<start1:
            #displace mid2 by start2 gaps at start2 to make pseudoknots match.
            for i in range(start1-start2):
                dict_seq2=insert_gap_seq_dict(dict_seq2, start2+mid_g2+1)
                  
            mid_d2=mid_d2+(start1-start2)
    
    trans_amnt=abs(mid_g1-mid_g2)

    diff1=mid_d1-mid_g1
    diff2=mid_d2-mid_g2
    
    if diff1==diff2:
        #only translation operation.
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

    if direction=="Left":
        #inversing the dictionnaries when aligning left hand side of pattern.
        dict_tba1=inverser_dico(dict_tba1)
        dict_tba2=inverser_dico(dict_tba2)
    
    if Pdiff > 0:
        # more pairings in 1.
        for i, elt2 in enumerate(cluster_2):
            if elt2 != cluster_1[i]:
                local_diff = elt2-cluster_1[i]
                if local_diff > 0:
                    # place local_diff gaps in 1 at pos cluster_1[i].
                    for k in range(local_diff):
                        dict_tba1 = insert_gap_seq_dict(dict_tba1, cluster_1[i])

                    # update cluster_1
                    for j in range(i, len(cluster_1)):
                        cluster_1[j] += local_diff
                else:
                    # place abs(local_diff) gaps in 2 at pos elt2.
                    for k in range(abs(local_diff)):
                        dict_tba2 = insert_gap_seq_dict(dict_tba2, elt2)
                    # update cluster_2
                    for j in range(i, len(cluster_2)):
                        cluster_2[j] += abs(local_diff)

    else:
        # more pairings in 2 or equal amount.
        for i, elt1 in enumerate(cluster_1):
            if elt1 != cluster_2[i]:
                local_diff = cluster_2[i]-elt1
                if local_diff > 0:
                    # place local_diff gaps in 1 at pos elt1.
                    for k in range(local_diff):
                        dict_tba1 = insert_gap_seq_dict(dict_tba1, elt1)
                    # update cluster_1
                    for j in range(i, len(cluster_1)):
                        cluster_1[j] += local_diff
                else:
                    # place abs(local_diff) gaps in 2 at pos cluster_2[i].
                    for k in range(abs(local_diff)):
                        dict_tba2 = insert_gap_seq_dict(
                            dict_tba2, cluster_2[i])
                    # update cluster_2
                    for j in range(i, len(cluster_2)):
                        cluster_2[j] += abs(local_diff)
                        
    if direction=="Left":      
        #reinversing the dictionnaries when finished aligning.
        dict_tba1=inverser_dico(dict_tba1)
        dict_tba2=inverser_dico(dict_tba2)

    #account for the rest of pairings by placing gaps.
    
    main_diff = len(dict_tba1)-len(dict_tba2)

    finish1=list(dict_tba1.keys())[-1]
    finish2=list(dict_tba2.keys())[-1]
    start1=list(dict_tba1.keys())[0]
    start2=list(dict_tba2.keys())[0]
    
    if main_diff>0:
        overhang_seq2=abs(main_diff)
        overhang_seq1=0
    elif main_diff<0:
        overhang_seq1=abs(main_diff)
        overhang_seq2=0
    else:
        overhang_seq2=0
        overhang_seq1=0
    
    if main_diff==0:
        return dict_tba1, dict_tba2, overhang_seq1, overhang_seq2
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
    
    return dict_tba1, dict_tba2, overhang_seq1, overhang_seq2

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
            
        dict1_L, dict2_L, left_overhang_pat1, left_overhang_pat2 = propagation_alignment(dict_tba1_L_translated, dict_tba2_L_translated, "Left")
    else:
        dict1_L, dict2_L, left_overhang_pat1, left_overhang_pat2 = propagation_alignment(dict_tba1_L, dict_tba2_L, "Left")
        
        start1=list(dict1_L.keys())[0]
        start2=list(dict2_L.keys())[0]
        
        dict1_L=dict_seq_translation(dict1_L,-start1)
        dict2_L=dict_seq_translation(dict2_L,-start2)
    
    # Now, translate the middle and right dictionnaries
    trans1 = list(dict1_L.values()).count("-") - start1
    trans2 = list(dict2_L.values()).count("-") - start2

    dict1_M = dict_seq_translation(dict_mid1, trans1)
    dict2_M = dict_seq_translation(dict_mid2, trans2)
    
    #aligning the right dictionnaries
    dict1_R, dict2_R, right_overhang_pat1, right_overhang_pat2 = propagation_alignment(dict_tba1_R, dict_tba2_R, "Right")

    #reagglomerating.
    seq1_L=dict_seq_reagglomerate(dict1_L)
    seq1_M=dict_seq_reagglomerate(dict1_M)
    seq1_R=dict_seq_reagglomerate(dict1_R)
    
    seq1=seq1_L+seq1_M+seq1_R
    
    seq2_L=dict_seq_reagglomerate(dict2_L)
    seq2_M=dict_seq_reagglomerate(dict2_M)
    seq2_R=dict_seq_reagglomerate(dict2_R)
    
    seq2=seq2_L+seq2_M+seq2_R

    return seq1, seq2, left_overhang_pat1, right_overhang_pat1, left_overhang_pat2, right_overhang_pat2

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
        seq1,seq2, left_overhang_pat1, right_overhang_pat1, left_overhang_pat2, right_overhang_pat2=inside_out_pat_alignment(pat1, pat2)
        
        pat1.left_overhang=left_overhang_pat1
        pat1.right_overhang=right_overhang_pat1
        
        pat2.left_overhang=left_overhang_pat2
        pat2.right_overhang=right_overhang_pat2
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
    Finds and marks with "O" and "C" the pairings that compose opening and closing overdivisions.
    """
    par="()"
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
        if dict_seq[elt]=="(":
            dict_seq[elt]="O"
        elif dict_seq[elt]==")":
            dict_seq[elt]="C"
            
    new_seq=''
    for elt in dict_seq.values():
        new_seq+=elt
    
    subdiv_list=[]
    subdiv_dict={}
    new_subdiv=True
    for i,elt in enumerate(new_seq):
        if elt=="O":
            if new_subdiv:
                new_subdiv=False
            subdiv_dict[i]="("

        elif elt =="C":
            if new_subdiv:
                new_subdiv=False
            subdiv_dict[i]=")"

        else:
            new_subdiv=True
            if subdiv_dict!={}:
                subdiv_list.append(subdiv_dict)
            subdiv_dict={}
            
    if subdiv_dict!={}:
        subdiv_list.append(subdiv_dict)
    
    return subdiv_list, new_seq

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

def surround(subdiv,sep):
    """
    Checks if an overdivision found is around a single pattern or not.
    """
    order_dict={}
    for i in range(len(subdiv)):
        order_dict[i]=None
    for i,subdiv_dict in enumerate(subdiv):
        for separator in sep:
            if separator.start <= list(subdiv_dict)[0] <= separator.finish and separator.start != -1:
                order_dict[i]=separator.nb

    order_before=order_dict[0]
    inter_order_dict=copy.deepcopy(order_dict)

    for nb,order in order_dict.items():
        if nb!=0:
            if order_before==order:
                del inter_order_dict[nb]
            else:
                order_before=order
    
    order_dict={}
    
    i=0
    for value in inter_order_dict.values():
        order_dict[i]=value
        i+=1
        
    if len(order_dict) <= 1:
        return False
    elif order_dict[0]==order_dict[1]-2 and len(order_dict)==2:
        return True
    else:
        return False

### CLASSES
                  
class Structure():
    """
    Class representing an aligned or not dotbracket sequence.
    
    It is formed of two lists : one with its patterns and one with its separators.
    """
    def __init__(self, sequence, ident=None, fam=None, AGU=None, subdiv_param=2):
        self.raw=sequence
        self.sequence=AGU
        self.length=len(sequence)
        self.subdiv_list, self.raw_nosubdiv=subdiv_finder(sequence, subdiv_param)
        count_par=sequence.count("(")+sequence.count(")")
        subdiv_count=0
        for elt in self.subdiv_list:
            subdiv_count+=len(elt)
        
        if subdiv_count == count_par:
            self.subdiv_list=[]
            self.raw_nosubdiv=self.raw
 
        sep,pat=slicer(self.raw_nosubdiv)

        if self.subdiv_list!=[]:
            if surround(self.subdiv_list,sep):
                self.raw_nosubdiv=self.raw
                self.subdiv_list=[]
                sep,pat=slicer(self.raw_nosubdiv)    
        
        self.separators=sep
        self.patterns=pat
        self.pattern_nb=len(pat)
        self.separator_nb=len(sep)
        self.isaligned=False
        self.id=ident
        self.family=fam
        self.alignedsequence=""
        self.alignedwith=None
    
    def __str__(self):
        tab=[["Type:","Structure"],
             ["Length:",self.length],
             ["ID:",self.id],
             ["Dotbracket Sequence:",self.raw],
             ["Subdiv Sequence:",self.raw_nosubdiv],
             ["Nucleotides Sequence:",self.sequence],
             ["Subdiv:",self.subdiv_list]]
        if self.isaligned:
            tab.append(["Aligned:","Yes"])
            tab.append(["Aligned Sequence:",self.alignedsequence])
            tab.append(["Aligned with:",self.alignedwith.alignedsequence])
        else:
            tab.append(["Aligned:","No"])
        tbp=tabulate(tab,tablefmt="fancy_outline",headers="firstrow")
        tbp+="\nSeparators: \n"
        for i,elt in enumerate(self.separators):
            tbp+="\n"+str(elt)
        tbp+="\nPatterns: \n"
        for i,elt in enumerate(self.patterns):
            tbp+="\n"+str(elt)
        return tbp

    def aligned(self):
        self.isaligned=True
        for pattern in self.patterns:
            if not pattern.isaligned:
                self.isaligned=False
        if not self.isaligned:
            print("Structure not yet aligned")
        
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
            for i,elt in enumerate(seqint):
                if elt == "O":
                    seq_dic[i]="("
                elif elt =="C":
                    seq_dic[i]=")"
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
        
        self.subdiv_index=raw.count('O')+raw.count('C')
        
        if self.subdiv_index!=0:
            self.subdiv_accounted=False
            
        self.pk_index_opened=raw.count("[")+raw.count("{")+raw.count("<")
        self.pk_index_closed=raw.count("]")+raw.count("}")+raw.count(">")
        self.pk_index=self.pk_index_closed+self.pk_index_opened
        
        if self.pk_index_closed!=0:
            self.pk_index_c_accounted=False
        if self.pk_index_opened!=0:
            self.pk_index_o_accounted=False

    def __str__(self):
        tab=[["Type:","Separator"],
             ["Length:",self.length],
             ["Order:",self.nb],
             ["Sequence:",self.sequence],
             ["Start:",self.start],
             ["Finish:",self.finish],
             ["Subdiv index:",self.subdiv_index]]
        return tabulate(tab,tablefmt="fancy_outline",headers="firstrow")

    
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
        self.paired=False
        self.subdiv_index=raw.count("O")+raw.count("C")
        self.pseudo_index=count_pseudo(raw)
        self.isaligned=False
        self.alignedwith=None
        self.alignedsequence=""
        self.left_overhang=0
        self.right_overhang=0
        
        
        self.pk_index_opened=raw.count("[")+raw.count("{")+raw.count("<")
        self.pk_index_closed=raw.count("]")+raw.count("}")+raw.count(">")
        self.pk_index=self.pk_index_closed+self.pk_index_opened
        
        if self.pk_index_closed!=0:
            self.pk_index_c_accounted=False
        
        if self.pk_index_opened!=0:
            self.pk_index_o_accounted=False
            
        
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
        
        tab=[["Type:","Pattern"],
             ["Length:",self.length],
             ["Order:",self.nb],
             ["Sequence:",self.sequence],
             ["Start:",self.start],
             ["Finish:",self.finish],
             ["Subdiv index:",self.subdiv_index],
             ["Left Overhang:",self.left_overhang],
             ["Right Overhang:",self.right_overhang]]
        if self.isaligned:
            tab.append(["Aligned:","Yes"])
            tab.append(["Aligned Sequence:",self.alignedsequence])
            tab.append(["Aligned with:",self.alignedwith.alignedsequence])
        else:
            tab.append(["Aligned:","No"])
        return tabulate(tab,tablefmt="fancy_outline",headers="firstrow")

EmptyPattern=Pattern('',[-1,-1],-1,-1)

### BASE SEPARATOR GAPS FUNCTIONS

def del_gaps(current_order, nb_gaps_deleted, order):
    """
    Updates the start and finish positions of further separators and patterns when gaps are deleted.
    """
    for elt in order:
        if elt.nb > current_order:
            if elt.start == -1:
                elt.finish-=nb_gaps_deleted
            else:
                elt.start-=nb_gaps_deleted
                elt.finish-=nb_gaps_deleted
    
def add_gaps(current_order,added_gaps,order_list):
    """
    Updates the start and finish positions of further separators and patterns when gaps are added.
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

### SEPARATOR ALIGNMENT FUNCTIONS

def pseudoknots_compensating(struct1, struct2, ordered1, ordered2, matching):
    """
    Spatially aware pseudoknots compensation inside of separators and non paired patterns.
    """
    #creating necessary separator and pattern matching 
    sep_opened_matching=[]
    sep_closed_matching=[]
    
    def gen_pair(triad, pairs):
        return triad[0][pairs[triad[1]].nb+triad[2]]
    

    def only_sep(upk, pairs):        
        if upk[0]=="o":
            if gen_pair(upk[1],pairs).pk_index_opened !=0 and gen_pair(upk[2],pairs).pk_index_opened !=0:
                if not(gen_pair(upk[1],pairs).pk_index_o_accounted) and not(gen_pair(upk[2],pairs).pk_index_o_accounted):
                    gen_pair(upk[1],pairs).pk_index_o_accounted = True
                    gen_pair(upk[2],pairs).pk_index_o_accounted = True
                    if upk[1][0]==ordered1:
                        return [gen_pair(upk[1],pairs),gen_pair(upk[2],pairs)]
                    else:
                        return [gen_pair(upk[2],pairs),gen_pair(upk[1],pairs)]
        
        elif upk[0]=="c":
            if gen_pair(upk[1],pairs).pk_index_closed !=0 and gen_pair(upk[2],pairs).pk_index_closed !=0:
                if not(gen_pair(upk[1],pairs).pk_index_c_accounted) and not(gen_pair(upk[2],pairs).pk_index_c_accounted):
                    gen_pair(upk[1],pairs).pk_index_c_accounted = True
                    gen_pair(upk[2],pairs).pk_index_c_accounted = True
                    if upk[1][0]==ordered1:
                        return [gen_pair(upk[1],pairs),gen_pair(upk[2],pairs)]
                    else:
                        return [gen_pair(upk[2],pairs),gen_pair(upk[1],pairs)]

        return False
            

    def non_paired(upk,pairs):
        if upk[4][0]=="b":
            bool_res = pairs[upk[1][0]].nb+upk[1][1] != 0
        elif upk[4][0]=="a":
            bool_res = pairs[upk[1][0]].nb+upk[1][1] != len(upk[4][1])-1
        if bool_res:
            if gen_pair(upk[2],pairs).alignedwith == EmptyPattern:
                if upk[0]=="o":
                    if gen_pair(upk[2],pairs).pk_index_opened != 0 and gen_pair(upk[3],pairs).pk_index_opened !=0:
                        if not(gen_pair(upk[3],pairs).pk_index_o_accounted) and not (gen_pair(upk[2],pairs).pk_index_o_accounted):
                            gen_pair(upk[3],pairs).pk_index_o_accounted = True
                            gen_pair(upk[2],pairs).pk_index_o_accounted = True
                            if upk[2][0]==ordered1:
                                return [gen_pair(upk[2],pairs),gen_pair(upk[3],pairs)]
                            else:
                                return [gen_pair(upk[3],pairs),gen_pair(upk[2],pairs)]
                elif upk[0]=="c":
                    if gen_pair(upk[2],pairs).pk_index_closed != 0 and gen_pair(upk[3],pairs).pk_index_closed !=0:
                        if not(gen_pair(upk[3],pairs).pk_index_c_accounted) and not (gen_pair(upk[2],pairs).pk_index_c_accounted):
                            gen_pair(upk[3],pairs).pk_index_c_accounted = True
                            gen_pair(upk[2],pairs).pk_index_c_accounted = True
                            return [gen_pair(upk[2],pairs),gen_pair(upk[3],pairs)]
                
        return False
    
    for pairs in matching:
        
        #BEFORE - OPENED    
        #test only separators BEFORE for opened pseudoknots
        
        unpacking=["o",[ordered1, 0, -1],[ordered2, 1, -1]]
        add=only_sep(unpacking,pairs)
        if add:
            sep_opened_matching.append(add)
    
        else:
            #test non paired Pattern BEFORE for opened pseudoknots in ordered1.
            unpacking=["o", [0,-1],[ordered1, 0, -2], [ordered2, 1, -1], ["b"]]
            add = non_paired(unpacking,pairs)
            if add:
                sep_opened_matching.append(add)
            
            #test non paired Pattern BEFORE for opened pseudoknots in ordered2.
            unpacking=["o", [1,-1],[ordered2, 1, -2], [ordered1, 0, -1], ["b"]]            
            add = non_paired(unpacking,pairs)
            if add:
                sep_opened_matching.append(add)
        
        #BEFORE - CLOSED
        #test only separators BEFORE for closed pseudoknots
        
        unpacking=["c",[ordered1, 0, -1],[ordered2, 1, -1]]
        add=only_sep(unpacking,pairs)
        if add:
            sep_closed_matching.append(add)
        
        else:
            #test non paired Pattern BEFORE for closed pseudoknots in ordered1
            unpacking=["c", [0,-1],[ordered1, 0, -2], [ordered2, 1, -1], ["b"]]
            add = non_paired(unpacking,pairs)
            if add:
                sep_closed_matching.append(add)
            
            unpacking=["c", [1,-1],[ordered2, 1, -2], [ordered1, 0, -1], ["b"]] 
            #test non paired Pattern BEFORE for closed pseudoknots in ordered2
            add = non_paired(unpacking,pairs)
            if add:
                sep_closed_matching.append(add)
        
        
        #AFTER - OPENED
        #test only separators AFTER for opened pseudoknots
        
        unpacking=["o",[ordered1, 0, 1],[ordered2, 1, 1]]
        add=only_sep(unpacking,pairs)
        if add:
            sep_opened_matching.append(add)
            
        else:
            #test non paired Pattern AFTER for opened pseudoknots in ordered1.
            unpacking=["o", [0,1],[ordered1, 0, 2], [ordered2, 1, 1], ["a",ordered1]]
            add=non_paired(unpacking,pairs)
            if add:
                sep_opened_matching.append(add)
            
            #test non paired Pattern AFTER for opened pseudoknots in ordered2.
            unpacking=["o", [1,1], [ordered2, 1, 2], [ordered1, 0, 1], ["a",ordered2]]
            add=non_paired(unpacking,pairs)
            if add:
                sep_opened_matching.append(add)
        
        #AFTER - CLOSED
        #test only separators AFTER for closed pseudoknots
        
        unpacking=["c",[ordered1, 0, 1],[ordered2, 1, 1]]
        add=only_sep(unpacking,pairs)
        if add:
            sep_closed_matching.append(add)
            
        else:
            #test non paired Pattern AFTER for closed pseudoknots in ordered1.
            unpacking=["c", [0,1],[ordered1, 0, 2], [ordered2, 1, 1], ["a",ordered1]]
            add = non_paired(unpacking,pairs)
            if add:
                sep_closed_matching.append(add)
            
            #test non paired Pattern AFTER for opened pseudoknots in ordered2.
            unpacking=["o", [1,1], [ordered2, 1, 2], [ordered1, 0, 1], ["a",ordered2]]        
            add = non_paired(unpacking,pairs)
            if add:
                sep_closed_matching.append(add)
    
    for pair in sep_opened_matching:
        dict_seq1={}
        for i,elt in enumerate(pair[0].sequence):
            dict_seq1[i]=elt
        
        dict_seq2={}
        for i,elt in enumerate(pair[1].sequence):
            dict_seq2[i]=elt
        
        start_diff1=0
        if pair[0].nb-1>=0 and isinstance(pair[0],Separator):
            if ordered1[pair[0].nb-1].alignedwith == EmptyPattern:
                start_diff1=ordered1[pair[0].nb-1].length
        
        start_diff2=0
        if pair[1].nb-1>=0 and isinstance(pair[1],Separator):
            if ordered2[pair[1].nb-1].alignedwith == EmptyPattern:
                start_diff2=ordered2[pair[1].nb-1].length
        
        opened="[{<"
        #determining the start of both pseudoknots
        start1=0
        while dict_seq1[start1] not in opened:
            start1+=1
            
        start2=0
        while dict_seq2[start2] not in opened:
            start2+=1
        
        start1+=start_diff1
        start2+=start_diff2
        
        #aligning if pk starts at different places in the separator
        if start1!=start2:
            
            diff=start1-start2
            
            if diff > 0:
                #start1>start2, # starts after in 1; placing gaps in 2.
                for i in range(abs(diff)):
                    dict_seq2 = insert_gap_seq_dict(dict_seq2,start2)
            elif diff <0:
                #start2>start1, # starts after in 2; placing gaps in 1.
                for i in range(abs(diff)):
                    dict_seq1 = insert_gap_seq_dict(dict_seq1,start1)
            
            if isinstance(pair[0],Pattern):
                pair[0].alignedsequence=dict_seq_reagglomerate(dict_seq1)
            else:
                pair[0].sequence=dict_seq_reagglomerate(dict_seq1)
                
            if isinstance(pair[1],Pattern):
                pair[1].alignedsequence=dict_seq_reagglomerate(dict_seq2)
            else:
                pair[1].sequence=dict_seq_reagglomerate(dict_seq2)
            
            nb_gaps1=pair[0].sequence.count('-')
            nb_gaps2=pair[1].sequence.count('-')
            
            pair[0].length+=nb_gaps1
            pair[1].length+=nb_gaps2
            
            pair[0].finish+=nb_gaps1
            pair[1].finish+=nb_gaps2

            add_gaps(pair[0].nb,nb_gaps1,ordered1)
            add_gaps(pair[1].nb,nb_gaps2,ordered2)

            struct1.length+=nb_gaps1
            struct2.length+=nb_gaps2
        
            
    for pair in sep_closed_matching:
    
        dict_seq1={}
        for i,elt in enumerate(pair[0].sequence):
            dict_seq1[i]=elt
        
        dict_seq2={}
        for i,elt in enumerate(pair[1].sequence):
            dict_seq2[i]=elt

        start_diff1=0
        if pair[0].nb-1>=0 and isinstance(pair[0],Separator):
            if ordered1[pair[0].nb-1].alignedwith == EmptyPattern:
                start_diff1=ordered1[pair[0].nb-1].length
        
        start_diff2=0
        if pair[1].nb-1>=0 and isinstance(pair[1],Separator):
            if ordered2[pair[1].nb-1].alignedwith == EmptyPattern:
                start_diff2=ordered2[pair[1].nb-1].length
        
        closed="]}>"
        #determining the start of both pseudoknots
        start1=0
        while dict_seq1[start1] not in closed:
            start1+=1
            
        start2=0
        while dict_seq2[start2] not in closed:
            start2+=1
        
        start1+=start_diff1
        start2+=start_diff2
        
        #aligning if pk starts at different places in the separator
        if start1!=start2:
            
            diff=start1-start2
            
            if diff > 0:
                #start1>start2, # starts after in 1; placing gaps in 2.
                for i in range(abs(diff)):
                    dict_seq2 = insert_gap_seq_dict(dict_seq2,start2)
            elif diff <0:
                #start2>start1, # starts after in 2; placing gaps in 1.
                for i in range(abs(diff)):
                    dict_seq1 = insert_gap_seq_dict(dict_seq1,start1)
            
            if isinstance(pair[0],Pattern):
                pair[0].alignedsequence=dict_seq_reagglomerate(dict_seq1)
            else:
                pair[0].sequence=dict_seq_reagglomerate(dict_seq1)
                
            if isinstance(pair[1],Pattern):
                pair[1].alignedsequence=dict_seq_reagglomerate(dict_seq2)
            else:
                pair[1].sequence=dict_seq_reagglomerate(dict_seq2)
                
            nb_gaps1=pair[0].sequence.count('-')
            nb_gaps2=pair[1].sequence.count('-')
            
            pair[0].length+=nb_gaps1
            pair[1].length+=nb_gaps2
            
            pair[0].finish+=nb_gaps1
            pair[1].finish+=nb_gaps2
            
            add_gaps(pair[0].nb,nb_gaps1,ordered1)
            add_gaps(pair[1].nb,nb_gaps2,ordered2)
            
            
            struct1.length+=nb_gaps1
            struct2.length+=nb_gaps2

def overdivision_compensating(struct1, struct2, ordered1, ordered2, matching):
    """
    Accounting for overdivision and aligning separators where overdivisions are present.
    
    We are only looking at separators present right before or after matched patterns.
    
    If there are spatial overdivision characters in both separators compared, we align them with eachother.
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

        #determining the start of both overdiv
        start1=0
        while dict_seq1[start1]!="O" and dict_seq1[start1]!="C":
            start1+=1
            
        start2=0
        while dict_seq2[start2]!="O" and dict_seq2[start2]!="C":
            start2+=1
        
        
        #aligning if subdiv starts at different places in the separator
        if start1!=start2:
            
            diff=start1-start2
            
            if diff > 0:
                #start1>start2, # starts after in sep1; placing gaps in sep2.
                for i in range(abs(diff)):
                    dict_seq2 = insert_gap_seq_dict(dict_seq2,start2)
                nb_gaps1=0
                nb_gaps2=abs(diff)
            elif diff <0:
                #start2>start1, # starts after in sep2; placing gaps in sep1.
                for i in range(abs(diff)):
                    dict_seq1 = insert_gap_seq_dict(dict_seq1,start1)
                nb_gaps1=abs(diff)
                nb_gaps2=0
            
            sep_pairs[0].sequence=dict_seq_reagglomerate(dict_seq1)
            sep_pairs[1].sequence=dict_seq_reagglomerate(dict_seq2)

            sep_pairs[0].length+=nb_gaps1
            sep_pairs[1].length+=nb_gaps2
            
            sep_pairs[0].finish+=nb_gaps1
            sep_pairs[1].finish+=nb_gaps2

            add_gaps(sep_pairs[0].nb,nb_gaps1,ordered1)
            add_gaps(sep_pairs[1].nb,nb_gaps2,ordered2)

            struct1.length+=nb_gaps1
            struct2.length+=nb_gaps2

def diff_overhang_calculator(pat1, pat2, struct1, struct2, order1, order2):
    """
    Function used to calculate the diff found in sep_gap_inserter.
    
    Takes into account the overhang of matched patterns and the length of separators.
    """
    
    if pat1.left_overhang !=0:
        # check for overhang deletion possibility on the left of pat1.
        if pat1.start > pat2.start:
            diff=pat1.start-pat2.start
            seq=pat1.alignedsequence
            if diff <= pat1.left_overhang:

                for i in range(diff):
                    seq=del_str(seq, 0)
                ret_diff=0
                    
                pat1.alignedsequence=seq
                    
                pat1.length-=diff
                pat1.finish-=diff
                struct1.length-=diff
                    
                del_gaps(pat1.nb, diff, order1)
            else:

                for i in range(pat1.left_overhang):
                    seq=del_str(seq, 0) 
                ret_diff = - (pat1.left_overhang - diff)
                    
                pat1.alignedsequence=seq
                    
                pat1.length-=pat1.left_overhang
                pat1.finish-=pat1.left_overhang
                struct1.length-=pat1.left_overhang
                    
                del_gaps(pat1.nb, pat1.left_overhang, order1)
        else:
            ret_diff=pat1.start-pat2.start
            
    elif pat2.left_overhang != 0:
        # check for overhang deletion possibility on the left of pat2
        if pat2.start > pat1.start:
            diff=pat2.start-pat1.start
            seq=pat2.alignedsequence
            
            if diff <= pat2.left_overhang:
                for i in range(diff):
                    seq=del_str(seq,0)

                ret_diff = 0
                
                pat2.alignedsequence=seq
                    
                pat2.length-=diff
                pat2.finish-=diff
                struct2.length-=diff
                    
                del_gaps(pat2.nb, diff, order2)
            else:

                for i in range(pat2.left_overhang):
                    seq=del_str(seq,0)
                ret_diff = pat2.left_overhang - diff
                
                pat2.alignedsequence=seq
                    
                pat2.length-=pat2.left_overhang
                pat2.finish-=pat2.left_overhang
                struct2.length-=pat2.left_overhang
                    
                del_gaps(pat2.nb, pat2.left_overhang, order2)
        else:
            ret_diff=pat1.start-pat2.start
    else:
        ret_diff=pat1.start-pat2.start  
            
    #now, check for the right sides:
    # if right overhang, check for (in this order) :
        # The end of the structure if no pk/od have already been accounted for.
        # A non matched pattern after.
        # A big separator after.
     
    if pat1.right_overhang !=0:
        seq=pat1.alignedsequence
        if order1[pat1.nb+1].nb==len(order1)-1:
            #end of the structure! completely remove overhang if no pk/od have been accounted for after!
            if order1[pat1.nb+1].subdiv_index!=0 and not(order1[pat1.nb+1].subdiv_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat1.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat1.alignedsequence=seq
                         
                pat1.length-=pat1.right_overhang
                pat1.finish-=pat1.right_overhang
                struct1.length-=pat1.right_overhang
                         
                del_gaps(pat1.nb, pat1.right_overhang, order1)
                
            elif (order1[pat1.nb+1].pk_index_closed !=0) and not(order1[pat1.nb+1].pk_index_c_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat1.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                
                pat1.alignedsequence=seq
                     
                pat1.length-=pat1.right_overhang
                pat1.finish-=pat1.right_overhang
                struct1.length-=pat1.right_overhang
                     
                del_gaps(pat1.nb, pat1.right_overhang, order1)
                
            elif (order1[pat1.nb+1].pk_index_opened!=0) and not (order1[pat1.nb+1].pk_index_o_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat1.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                
                pat1.alignedsequence=seq
                     
                pat1.length-=pat1.right_overhang
                pat1.finish-=pat1.right_overhang
                struct1.length-=pat1.right_overhang
                     
                del_gaps(pat1.nb, pat1.right_overhang, order1)
             
        elif order1[pat1.nb+2].alignedwith == EmptyPattern:
            #take the length of the empty pattern and reduce the gaps accordingly.
            full_length=(order1[pat1.nb+2].length + order1[pat1.nb+1].length + order1[pat1.nb+3].length)
            diff = pat1.right_overhang - full_length
            
            if diff >= 0:
                #overhang > length of pattern + length of separators, reduce by length of all. (very unlikely)
                back = len(seq)-1
                cpt=0
                while cpt <= full_length-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat1.alignedsequence=seq

                pat1.length-=full_length
                pat1.finish-=full_length
                struct1.length-=full_length
                     
                del_gaps(pat1.nb, full_length, order1)
            else:
                #reduce by the amount of overhang.
                back = len(seq)-1
                cpt=0
                while cpt <= pat1.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat1.alignedsequence=seq
                     
                pat1.length-=pat1.right_overhang
                pat1.finish-=pat1.right_overhang
                struct1.length-=pat1.right_overhang
                     
                del_gaps(pat1.nb, pat1.right_overhang, order1)
        else:
            #separator after!
            diff = pat1.right_overhang - order1[pat1.nb+1].length
            if diff >=0:
                #overhang bigger than separator! Reduce by length of separator.
                back = len(seq)-1
                cpt=0
                if order1[pat1.nb+1].length != 0:
                    while cpt <= order1[pat1.nb+1].length -1 :
                        seq=del_str(seq,back)
                        cpt+=1
                        back-=1
    
                        
                    pat1.alignedsequence=seq
                         
                    pat1.length-=order1[pat1.nb+1].length
                    pat1.finish-=order1[pat1.nb+1].length
                    struct1.length-=order1[pat1.nb+1].length
                         
                    del_gaps(pat1.nb, order1[pat1.nb+1].length, order1)
            else:
                #seprator bigger!, reduce by length of overhang!
                back = len(seq)-1
                cpt=0
                while cpt <= pat1.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat1.alignedsequence=seq
                     
                pat1.length-=pat1.right_overhang
                pat1.finish-=pat1.right_overhang
                struct1.length-=pat1.right_overhang
                     
                del_gaps(pat1.nb, pat1.right_overhang, order1)
                
    elif pat2.right_overhang !=0:
        seq=pat2.alignedsequence
        if order2[pat2.nb+1].nb==len(order2)-1:
            #end of the structure! completely remove overhang if no pk/od have been accounted for after!
            if order2[pat2.nb+1].subdiv_index!=0 and not(order2[pat2.nb+1].subdiv_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat2.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                
                pat2.alignedsequence=seq
                     
                pat2.length-=pat2.right_overhang
                pat2.finish-=pat2.right_overhang
                struct2.length-=pat2.right_overhang
                     
                del_gaps(pat2.nb, pat2.right_overhang, order2)
                
            elif (order2[pat2.nb+1].pk_index_closed !=0) and not(order2[pat2.nb+1].pk_index_c_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat2.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                
                pat2.alignedsequence=seq
                     
                pat2.length-=pat2.right_overhang
                pat2.finish-=pat2.right_overhang
                struct2.length-=pat2.right_overhang
                     
                del_gaps(pat2.nb, pat2.right_overhang, order2)
                
            elif (order2[pat2.nb+1].pk_index_opened!=0) and not (order2[pat2.nb+1].pk_index_o_accounted):
                back = len(seq)-1
                cpt=0
                while cpt <= pat2.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                
                pat2.alignedsequence=seq
                     
                pat2.length-=pat2.right_overhang
                pat2.finish-=pat2.right_overhang
                struct2.length-=pat2.right_overhang
                     
                del_gaps(pat2.nb, pat2.right_overhang, order2)
                
        elif order2[pat2.nb+2].alignedwith == EmptyPattern:
            #take the length of the empty pattern and reduce the gaps accordingly.
            full_length=(order2[pat2.nb+2].length + order2[pat2.nb+1].length + order2[pat2.nb+3].length)
            diff = pat2.right_overhang - full_length
           
            if diff >= 0:
                #overhang > length of pattern + length of separators, reduce by length of all. (very unlikely)
                back = len(seq)-1
                cpt=0
                while cpt <= full_length-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat2.alignedsequence=seq
                     
                pat2.length-=full_length
                pat2.finish-=full_length
                struct2.length-=full_length
                     
                del_gaps(pat2.nb, full_length, order2)
            else:
                #reduce by the amount of overhang.
                back = len(seq)-1
                cpt=0
                while cpt <= pat2.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat2.alignedsequence=seq
                     
                pat2.length-=pat2.right_overhang
                pat2.finish-=pat2.right_overhang
                struct2.length-=pat2.right_overhang
                     
                del_gaps(pat2.nb, pat2.right_overhang, order2)
        else:
            #separator after!
            diff = pat2.right_overhang - order2[pat2.nb+1].length
            
            if diff >=0:
                #overhang bigger than separator! Reduce by length of seprator.
                back = len(seq)-1
                if order2[pat2.nb+1].length != 0:
                    cpt=0
                    while cpt <= order2[pat2.nb+1].length-1:
                        seq=del_str(seq,back)
                        cpt+=1
                        back-=1
                        
                    pat2.alignedsequence=seq
                         
                    pat2.length-=order2[pat2.nb+1].length
                    pat2.finish-=order2[pat2.nb+1].length
                    struct2.length-=order2[pat2.nb+1].length
                         
                    del_gaps(pat2.nb, order2[pat2.nb+1].length, order2)
            else:
                #seprator bigger!, reduce by length of overhang!
                back = len(seq)-1
                cpt=0
                while cpt <= pat2.right_overhang-1:
                    seq=del_str(seq,back)
                    cpt+=1
                    back-=1
                    
                pat2.alignedsequence=seq
                     
                pat2.length-=pat2.right_overhang
                pat2.finish-=pat2.right_overhang
                struct2.length-=pat2.right_overhang
                     
                del_gaps(pat2.nb, pat2.right_overhang, order2)
                
    return ret_diff

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
            
            diff = diff_overhang_calculator(elt[0], elt[1], struct1, struct2, ordered1, ordered2)
            
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

def discriminate_po_priority(struct1, struct2, matching):
    """
    This function checks whether this situation is present inside separators:
        
    #####[[[[[
    [[[[[#####
         
         OR
    
    [[[[[#####
    #####[[[[[

    with #### the overdivision and [[[[ the pseudoknots.
    
    In that case, the entirety of the bloc should be considered for alignment and ignored by both function after.
    
    Returns a boolean.
    In all cases, pseudoknots_compensating should be executed since it tests for pseudoknots also in non paired patterns.
    ALL OVERDIVISION SHOULD BE IN SEPARATORS.
    """
    
    pk=["[","{","<",">","}","]"]
    # first, create the matched separators list:
    order1=struct1.order_list()
    order2=struct2.order_list()
    
    matched_sep=[]
    accounted_sep=[]

    for tup_pat in matching:
        pat1=tup_pat[0]
        pat2=tup_pat[1]
        #Before the two matched patterns:
        if (order1[pat1.nb-1],order2[pat2.nb-1]) not in accounted_sep:
            matched_sep.append([order1[pat1.nb-1],order2[pat2.nb-1]])
            accounted_sep.append((order1[pat1.nb-1],order2[pat2.nb-1]))
        #After the two matched patterns:
        if (order1[pat1.nb+1],order2[pat2.nb+1]) not in accounted_sep:
            matched_sep.append([order1[pat1.nb+1],order2[pat2.nb+1]])
            accounted_sep.append((order1[pat1.nb+1],order2[pat2.nb+1]))
         
    criss_cross = []
    #Now testing for overdivision in known matched separators:
    for sepmatch in matched_sep:
        sep1=sepmatch[0]
        sep2=sepmatch[1]
        if sep1.subdiv_index!=0 and sep2.subdiv_index!=0:
            #overdivision are present.
            if sep1.pk_index!=0 and sep2.pk_index!=0:
                #pseudoknots are also present.
                #check if pseudoknots and overdivision have a criss-cross pattern and mark the separators for alignment.
                if sep1.pk_index_closed!=0:
                    start_od1=np.char.find(sep1.sequence, "C")
                else:
                    start_od1=np.char.find(sep1.sequence, "O")
                
                if sep2.pk_index_closed!=0:
                    start_od2=np.char.find(sep2.sequence, "C")
                else:
                    start_od2=np.char.find(sep2.sequence, "O")
                    
                pk_array1=[]
                pk_array2=[]
                for elt in pk:
                    pk_array1.append(np.char.find(sep1.sequence,elt))
                    pk_array2.append(np.char.find(sep2.sequence,elt))
                
                for i,elt in enumerate(pk_array1):
                    if elt == -1:
                        pk_array1[i]=np.inf
                
                for i,elt in enumerate(pk_array2):
                    if elt == -1:
                        pk_array2[i]=np.inf
                
                start_pk1=min(pk_array1)
                start_pk2=min(pk_array2)

                if start_pk1 < start_od1 and start_pk2 > start_od2:
                    #criss-cross detected! - mark the considered separators!
                    bk_start1 = start_pk1
                    bk_start2 = start_od2
                    criss_cross.append([sep1,sep2, bk_start1, bk_start2])
                elif start_pk1 > start_od1 and start_pk2 < start_od2:
                    #criss-cross detected! - mark the considered separators!
                    bk_start1 = start_od1
                    bk_start2 = start_pk2
                    criss_cross.append([sep1,sep2, bk_start1, bk_start2])
           
    #Aligning the start of each blocks
    for sep_mark in criss_cross:
        sep1=sep_mark[0]
        sep2=sep_mark[1]
        start1=int(sep_mark[2])
        start2=int(sep_mark[3])
        
        
        dict_seq1={}
        for i,elt in enumerate(sep1.sequence):
            dict_seq1[i]=elt
        
        dict_seq2={}
        for i,elt in enumerate(sep2.sequence):
            dict_seq2[i]=elt
         
        if sep1.pk_index_closed!=0:
            sep1.pk_index_c_accounted=True
        if sep1.pk_index_opened!=0:
            sep1.pk_index_o_accounted=True
        sep1.subdiv_accounted=True
            
        if sep2.pk_index_closed!=0:
            sep2.pk_index_c_accounted=True
        if sep2.pk_index_opened!=0:
            sep2.pk_index_o_accounted=True
        sep2.subdiv_accounted=True
            
        if start1!=start2:
            diff=start1-start2
            
            if diff > 0:
                #start1>start2, # starts after in sep1; placing gaps in sep2.
                for i in range(abs(diff)):
                    dict_seq2 = insert_gap_seq_dict(dict_seq2,start2)
                nb_gaps1=0
                nb_gaps2=abs(diff)
            elif diff <0:
                #start2>start1, # starts after in sep2; placing gaps in sep1.
                for i in range(abs(diff)):
                    dict_seq1 = insert_gap_seq_dict(dict_seq1,start1)
                nb_gaps1=abs(diff)
                nb_gaps2=0
            
            sep1.sequence=dict_seq_reagglomerate(dict_seq1)
            sep2.sequence=dict_seq_reagglomerate(dict_seq2)

            sep1.length+=nb_gaps1
            sep2.length+=nb_gaps2
            
            sep1.finish+=nb_gaps1
            sep2.finish+=nb_gaps2

            add_gaps(sep1.nb,nb_gaps1,order1)
            add_gaps(sep2.nb,nb_gaps2,order2)

            struct1.length+=nb_gaps1
            struct2.length+=nb_gaps2
           
### MATCHING FUNCTION

def pair_pat_score_pass1(pat1,pat2):
    """
    Gives a pairing score for first pass pattern recognition based on the aptamat distance, the number of pairings and the length difference.
    """
    apta_dist=AF.compute_distance_clustering(AF.SecondaryStructure(pat1.sequence),AF.SecondaryStructure(pat2.sequence), "cityblock", "slow")
    
    pat1_sc = (pat1.sequence.count("(")+pat1.sequence.count(")"))
    pat2_sc = (pat2.sequence.count("(")+pat2.sequence.count(")"))
    
    #pseudo=abs(pat1.pseudo_index - pat2.pseudo_index)
    
    return abs(pat1_sc - pat2_sc) + apta_dist + abs(pat1.length - pat2.length)# + pseudo

def pair_pat_score_pass2(pat1,pat2):
    """
    Gives a pairing score for second pass pattern recognition based on the aptamat distance and the number of pairings.
    """
    apta_dist=AF.compute_distance_clustering(AF.SecondaryStructure(pat1.sequence),AF.SecondaryStructure(pat2.sequence), "cityblock", "slow")
    
    pat1_sc = (pat1.sequence.count("(")+pat1.sequence.count(")"))
    pat2_sc = (pat2.sequence.count("(")+pat2.sequence.count(")"))
    
    #pseudo=abs(pat1.pseudo_index - pat2.pseudo_index)

    return abs(pat1_sc - pat2_sc) + apta_dist # + pseudo

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
            
            for elt in matching_pass3:
                elt[0].paired=True
                elt[1].paired=True
            for pat1 in struct1.patterns:
                if not pat1.paired:
                    pat1.aligned(None,pat1.sequence)

            for pat2 in struct2.patterns:
                if not pat2.paired:
                    pat2.aligned(None,pat2.sequence)

            return matching_pass3
        else:
            if verbose:
                print("Pass2 valid, returning matching.")

            for elt in matching_pass2:
                elt[0].paired=True
                elt[1].paired=True
            
            for pat1 in struct1.patterns:
                if not pat1.paired:
                    pat1.aligned(None,pat1.sequence)
            
            for pat2 in struct2.patterns:
                if not pat2.paired:
                    pat2.aligned(None,pat2.sequence)
                    
            return matching_pass2
    else:
        if verbose:
            print("Pass1 valid, returning matching.")

        for elt in matching_pass1:
            elt[0].paired=True
            elt[1].paired=True
        
        for pat1 in struct1.patterns:
            if not pat1.paired:
                pat1.aligned(None,pat1.sequence)
                
        for pat2 in struct2.patterns:
            if not pat2.paired:
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
        
        - Testing for compatibility using a pairing score and three passes.
        
        - If the best match is not compatible with an alignment, use naive matching instead.
        
        ALIGNING
        
        - Using inside-out alignment to align each pattern with its match.
            - mark them as aligned and input each other in the `self.alignedwith` variable.
        
        - Test for overdivisions alignment possibilities and aligns the separators compatible.
        
        - When done, match the length of the two structures with gaps placed in the separators where it minimizes the distance.
        
        - Mark the structures as aligned and input each other in the `self.alignedwith` variable.
        
        - return the structures and the new distance.
    ____________________
    
    """

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
            print("Same structures.")
    else:
        if verbose:
            print("\nDetermining the optimal pattern match if possible \n")
        try:
            matching = matching_finder(struct1, struct2, verbose)
        except Exception as e:
            e = f"An error occured when matching patterns on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
            raise MatchingError(e)
    
    if verbose:
        print("\nAligning patterns.")
    

    for elt in matching:
        order1=struct1.order_list()
        order2=struct2.order_list()
        try:
            pattern_alignment(struct1,struct2,elt[0],elt[1],order1,order2,verbose)
        except Exception as e:
            e = f"An error occured when aligning patterns on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
            raise PatternAlignmentError(e)
        
    if verbose:
        print("\nAccounting for pseudoknots and overdivision.")

    order1=struct1.order_list()
    order2=struct2.order_list()
    
    try:
        discriminate_po_priority(struct1, struct2, matching)
    except Exception as e:
        e=f"An error occured when checking for overdiv/pseudoknots priority on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
        raise OPKPrioError(e)
        
    try:
        pseudoknots_compensating(struct1, struct2, order1, order2, matching)    
    except Exception as e:
        e=f"An error occured when compensating for pseudoknots on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
        raise PKCompensatingError(e)
            
    try:
        overdivision_compensating(struct1, struct2, order1, order2, matching)
    except Exception as e:
        e = f"An error occured when compensating for overdivision on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
        raise ODCompensatingError(e)
        
    if verbose:
        print("\nAdding gaps in separators for length and pattern matching")

    try:
        separator_compensating(struct1, struct2, matching)
    except Exception as e:
        e = f"An error occured when aligning separators on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
        raise SepCompensatingError(e)
    
    if verbose:
        print("\nUpdating last parameters and finishing\n")
    struct1.alignedwith=struct2
    struct2.alignedwith=struct1
    
    struct1.aligned()
    struct2.aligned()
  
    struct1.alignedsequence = struct1.reagglomerate()
    struct2.alignedsequence = struct2.reagglomerate()
    
    return struct1,struct2

#_____________________________MAIN FUNCTIONS_____________________________

def opt_subdiv(seq1, seq2, depth):
    """
    opt_subdiv finds the best subdiv parameter via brute force. This function runs the alignments in parallel.
    """
    if mp.cpu_count() < depth:
        try:
            pool = mp.Pool(mp.cpu_count())
        except Exception as e:
            e = f"Error creating multiprocessing pool: \n     {e}"
            raise MultiprocessingError(e)
    else:
        try:
            pool = mp.Pool(depth)
        except Exception as e:
            e = f"Error creating multiprocessing pool: \n     {e}"
            raise MultiprocessingError(e)
        
    structure_list=[]
    for i in range(1,depth+1):
        try:
            structure_list.append([Structure(seq1,subdiv_param=i),Structure(seq2,subdiv_param=i)])
        except Exception as e:
            e = f"Error when creating structures on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
            raise StructBuildError(e)
    
    results=[]
    try:
        for struct1, struct2 in pool.starmap(full_alignment, [(structs[0], structs[1]) for structs in structure_list]):
            results.append([struct1,struct2])
    except Exception as e:
        e = f"An error occured while aligning in parallel: \n     {e}"
        pool.close()
        raise GeneralAlignmentError(e)
    
    dists=[]
    
    
    for struct1, struct2 in results:
        dists.append(AF.compute_distance_clustering(AF.SecondaryStructure(struct1.alignedsequence),AF.SecondaryStructure(struct2.alignedsequence), "cityblock", "slow"))
    
    min_dist = np.argmin(dists)
    
    return results[min_dist][0], results[min_dist][1]

def clustering_opt_subdiv(seq1,seq2, depth,ident1=None, fam1=None, AGU1=None,ident2=None, fam2=None, AGU2=None):
    """
    clustering_opt_subdiv finds the best subdiv parameter via brute force. This function does NOT run the alignments in parallel.
    """
    structure_list=[]
    
    try:
        for i in range(1,depth+1):
            structure_list.append([Structure(seq1,subdiv_param=i,ident=ident1,fam=fam1,AGU=AGU1),Structure(seq2,subdiv_param=i,ident=ident2,fam=fam2,AGU=AGU2)])
    except Exception as e:
        e = f"Error when creating Structures on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
        raise StructBuildError(e)
        
    results=[]
    for struct1, struct2 in structure_list:
        try:
            results.append(full_alignment(struct1,struct2))
        except Exception as e:
            e = f"An error occured while aligning on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}"
            raise GeneralAlignmentError(e)
    dists=[]
    
    for struct1, struct2 in results:
        dists.append(AF.compute_distance_clustering(AF.SecondaryStructure(struct1.alignedsequence),AF.SecondaryStructure(struct2.alignedsequence), "cityblock", "slow"))
    
    min_dist = np.argmin(dists)
    
    return results[min_dist][0], results[min_dist][1]
    
def main():

    parser = argparse.ArgumentParser(description="AptAlign is an alignment algorithm designed around pattern recognition.\n"
                                                 "Use -s to input only two structures directly in the command line.\n"
                                                 "Use -v to toggle verbose mode.\n"
                                                 "Use -d to set the depth value for optimal overdivision parameter search.\n"
                                                 "Use -u to toggle the unoptimised version of the optimal overdivision parameter search.")

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
    
    parser.add_argument('-d',
                        '--depth',
                        nargs=1,
                        type=int,
                        default=[10],
                        help='Search depth for optimal overdivision parameter.')
                        
    parser.add_argument('-u',
                        '--unoptimised',
                        help="Use the unoptimised overdivision parameter calculator function.",
                        default=False,
                        action="store_true")
    try:
        args = parser.parse_args()
        
    except Exception as e:
        print(f"Error parsing arguments: {e}")
        sys.exit(1)
        
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if args.structures is not None:
        a=time.time()
        
        if args.unoptimised:
            struct1, struct2 = clustering_opt_subdiv(args.structures[0],args.structures[1],args.depth[0])
        else:
            struct1, struct2 = opt_subdiv(args.structures[0],args.structures[1],args.depth[0])
        
        initial_dist=AF.compute_distance_clustering(AF.SecondaryStructure(struct1.raw),AF.SecondaryStructure(struct2.raw), "cityblock", "slow")
        new_dist = AF.compute_distance_clustering(AF.SecondaryStructure(struct1.alignedsequence),AF.SecondaryStructure(struct2.alignedsequence), "cityblock", "slow")
        
        str1=struct1.alignedsequence
        str2=struct2.alignedsequence
        
        if not args.verbose:
            
            if initial_dist!=0:
                improvement = round((initial_dist-new_dist)/initial_dist *100,2)
            else:
                improvement = 1
            b=time.time()
            
            for i in range(int(len(str1)/190)+1):
                str1=insert_str(str1,i*190, "\n")
            for i in range(int(len(str2)/190)+1):
                str2=insert_str(str2,i*190, "\n")
                count=i
            
            tab=[["Structure 1:"+(count+2)*"\n"+"Structure 2:",str1+"\n"+str2],
                 ["Improvement: ",str(initial_dist)+" -> "+str(new_dist)+" | in %: "+str(improvement)+"%"]]
            print("\n")
            print(tabulate(tab,headers=["Results","In "+str(b-a)+"s"],tablefmt="fancy_grid"))
            print("\n")
        else:
            if initial_dist!=0:
               improvement = round((initial_dist-new_dist)/initial_dist *100,2)
            else:
               improvement = 1
            b=time.time()
            
            for i in range(int(len(str1)/190)+1):
                str1=insert_str(str1,i*190, "\n")
            for i in range(int(len(str2)/190)+1):
                str2=insert_str(str2,i*190, "\n")
                count=i
            
            tab=[["Structure 1:"+(count+2)*"\n"+"Structure 2:",str1+"\n"+str2],
                 ["Improvement: ",str(initial_dist)+" -> "+str(new_dist)+" | in %: "+str(improvement)+"%"]]
            print("\n")
            print(tabulate(tab,headers=["Results","In "+str(b-a)+"s"],tablefmt="fancy_grid"))
            print("\n")
            print("\nDetailed results:\n")
            print("______________________________________\n")
            print("First Structure","\n")
            print(struct1)
            print("______________________________________\n")
            print("Second Structure","\n")
            print(struct2)
        sys.exit(0)
    else:
        raise ValueError("The Structures were incorrectly parsed from command line.")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"General Alignment Error on line {sys.exc_info()[2].tb_next.tb_lineno}: \n     {e}")
        sys.exit(1)
