import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptafast'))
sys.path.append(root_path)

#import numpy as np
import AptaFast as AF

import math as m
#import time
#import multiprocessing


### REUSING THE FUNCTION FROM `clustering_AptaMat.py` USING THE SAME FILE TYPE AND CONSTRUCTION

def initialize_dataset(structure_file):
    structure_list = []
    family = {}
    with open(structure_file, 'r') as file:
        for line in file:
            content = line.strip().split()
            if content:
                # print(line)
                if line.startswith('FAMILY'):
                    pass
                else:
                    try:
                        family[content[0]] += 1
                    except KeyError:
                        family[content[0]] = 1
    
                if AF.Dotbracket.is_dotbracket(content[3]):
                    structure = AF.SecondaryStructure(dotbracket=content[3], sequence=content[2],
                                                           id=content[1].split('.')[0])
                    # AptaMat._create_fasta(structure)
                    structure_list.append(structure)
    return structure_list,family



def double_mode_init(struct1,struct2):
    
    if AF.Dotbracket.is_dotbracket(struct1):
        structure_1 = AF.SecondaryStructure(dotbracket=struct1)
    if AF.Dotbracket.is_dotbracket(struct2):
        structure_2 = AF.SecondaryStructure(dotbracket=struct2)
    
    return structure_1,structure_2


    
def del_str(string,num):
    string=list(string)
    del(string[num])
    res=""
    for elt in string:
        res+=elt
    return res
    
def find_dash(string):
    res=[]
    string = list(string)
    for i,elt in enumerate(string):
        if elt =="-":
            res.append(i)
    return res
    
def insert_str(string,num,char):
    string=list(string)
    string.insert(num,char)
    res=""
    for elt in string:
        res+=elt
    return res

def double_mode_brute_force(struct1,struct2,max_size=0):
    
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)

    
    fill1=max_size-len(struct1)
    fill2=max_size-len(struct2)
    
    struct1_prealign=fill1*"-"+struct1
    struct2_prealign=fill2*"-"+struct2
    
    # on place le 1er gap -> len(struct1)+fill1 positions possibles ; on place le second gap -> len(struct1)+fill1-1 positions possibles : max_size factorielles/fill1 factorielle positions possibles.
    i1 = int(m.factorial(max_size)/m.factorial(fill1))
    i2 = int(m.factorial(max_size)/m.factorial(fill2))
    if find_dash(struct2_prealign) == []:
        i2=0
    if find_dash(struct1_prealign) == []:
        i1=0
    all_struct1=[struct1_prealign]
    all_struct2=[struct2_prealign]
    
    temp_struct1=struct1_prealign
    ins=0
    comb=0
    print(i1,i2)
    for i in range(i1):
        dashes=find_dash(temp_struct1)
        #print(dashes)
        temp_struct1=del_str(temp_struct1,dashes[0])
        temp_struct1=insert_str(temp_struct1,ins,"-")
        if temp_struct1 not in all_struct1: 
            all_struct1.append(temp_struct1)
        ins+=1
        if ins == len(temp_struct1):
            comb+=1
            if comb == len(temp_struct1):
                comb=0
            ins=0
    
    temp_struct2=struct2_prealign
    ins=0
    comb=0
    for i in range(i2):
        dashes=find_dash(temp_struct2)
        print(dashes)
        temp_struct2=del_str(temp_struct2,dashes[0])
        temp_struct2=insert_str(temp_struct2,ins,"-")
        
        if temp_struct2 not in all_struct2: 
            all_struct2.append(temp_struct2)
        ins+=1
        if ins == len(temp_struct2):
            comb+=1
            if comb == len(temp_struct2):
                comb=0
            ins=0
    
    return all_struct1, all_struct2

    
def brute_force_calc(l_struct1,l_struct2):
    dict_dist={}
    for elt1 in l_struct1:
        for elt2 in l_struct2:
            dict_dist[AF.compute_distance_clustering(AF.SecondaryStructure(elt1),AF.SecondaryStructure(elt2), "cityblock", "slow")]=(elt1,elt2)
    keep=min(dict_dist.keys())
    return dict_dist[keep],keep

