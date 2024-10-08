import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptafast'))
sys.path.append(root_path)

#import numpy as np
import AptaFast as AF
import tqdm
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


def arrangement(struct,max_size):
    
    fill=max_size-len(struct)
    
    struct_prealign=fill*"-"+struct
    
    #i=int(m.factorial(max_size)/m.factorial(fill))
    
    if find_dash(struct_prealign) == []:
        return [struct_prealign]
    
    all_struct=[struct_prealign]
    
    temp_struct=struct_prealign
    dashes=find_dash(temp_struct)
    for i in range(fill):
        for j,elt in enumerate(all_struct):
            dashes=find_dash(elt)
            for elt1 in one_range(dashes,i,elt,all_struct):
                all_struct.append(elt1)
    return all_struct

def one_range(dashes,num,curr_struct,all_struct):
    curr_dash=dashes[num]
    res=[]
    temp_struct=curr_struct
    temp_struct=del_str(temp_struct,curr_dash)
    temp_struct=insert_str(temp_struct,0,"-")
    if temp_struct not in all_struct:
        if temp_struct not in res:
            res.append(temp_struct)
    for i in range(1,len(curr_struct)):
        temp_struct=del_str(temp_struct,i-1)
        temp_struct=insert_str(temp_struct,i,"-")
        if temp_struct not in all_struct:
            if temp_struct not in res:
                res.append(temp_struct)
    return res

def brute_force_calc(struct1,struct2,max_size=0):
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)
    print("Generating arrangement")
    l_struct1=arrangement(struct1,max_size)
    l_struct2=arrangement(struct2,max_size)
    dict_dist={}
    print("Calculating")
    for elt1 in l_struct1:
        for elt2 in l_struct2:
            dict_dist[AF.compute_distance_clustering(AF.SecondaryStructure(elt1),AF.SecondaryStructure(elt2), "cityblock", "slow")]=(elt1,elt2)
    keep=min(dict_dist.keys())
    return dict_dist[keep],keep

