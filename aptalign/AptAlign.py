import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptafast'))
sys.path.append(root_path)

import AptaFast as AF
#import time
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np

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


def arrangement(struct,max_size):
    """This function calculates all the possible gap placement of one structure.
    
    *struct* is the base structure.
    *max_size* is the size of the alignment.
    
    """
    fill=max_size-len(struct)
    
    struct_prealign=fill*"-"+struct
    
    #i=int(m.factorial(max_size)/m.factorial(fill))
    
    if find_dash(struct_prealign) == []:
        return [struct_prealign]
    
    all_struct={struct_prealign:False}
    print("PASS1")
    temp_struct=struct_prealign
    dashes=find_dash(temp_struct)
    for i in range(fill):
        print(round((i/fill)*100,3),"%")
        for j,elt in enumerate(list(all_struct)):
            if not all_struct[elt]:
                dashes=find_dash(elt)
                all_struct[elt] = True
                res1=one_range(dashes,i,elt,all_struct)
                for elt1 in res1:
                    all_struct[elt1]=False
    
    struct_prealign=struct+fill*"-"
    
    if find_dash(struct_prealign) == []:
        return [struct_prealign]
    
    all_struct[struct_prealign]=False
    print("PASS2")
    temp_struct=struct_prealign
    dashes=find_dash(temp_struct)
    for i in range(fill):
        print(round((i/fill)*100,3),"%")
        for j,elt in enumerate(list(all_struct)):
            if not all_struct[elt]:
                dashes=find_dash(elt)
                all_struct[elt] = True
                res1=one_range(dashes,i,elt,all_struct)
                for elt1 in res1:
                    all_struct[elt1]=False
    
    return all_struct

def one_range(dashes,num,curr_struct,all_struct):
    """returns a 'swoop' of gapped structures with respect to an existing gap and a current structure
    
        *dashes* is the list of positions of all gaps from curr_struct.
        *num* is the gap number the swoop is made from
        *curr_struct* is the structure from which the swoop is calculated
        *all_struct* is the dictionnary of all already calculated structures.
        
        returns a list of new structures not already calculated
    """
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

def calc_dist(elt1,l_struct2):
    """Module used in the multiprocessing to calculate a batch of distances"""
    res=[]
    for elt2 in l_struct2:
        res.append((AF.compute_distance_clustering(AF.SecondaryStructure(elt1),AF.SecondaryStructure(elt2), "cityblock", "slow"),elt1,elt2))
    return res

def brute_force_calc(struct1,struct2,max_size=0):
    """Function used to align two structures using a burte force method
    
    This calculates all possible alignments and returns the one with the smallest AptaMat distance.
    
    max_size is the size the structure will be aligned to.
        Is equal to the length of the biggest structure if under.
    """
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)
    print("Generating arrangement")
    l_struct1=arrangement(struct1,max_size)
    l_struct2=arrangement(struct2,max_size)
    dict_dist={}
    print("S1: ",len(l_struct1),"| S2: ", len(l_struct2))
    print("This is a multiprocessed program, you have",mp.cpu_count(),"cores in your CPU.")
    nb=int(input("How much do you want to use? "))
    print("Creating pool on",nb,"cores.\n")
    print("Working...\n")
    pooling=mp.Pool(nb)
    
    # for elt1 in l_struct1:
    #     for elt2 in l_struct2:
    #         dict_dist[AF.compute_distance_clustering(AF.SecondaryStructure(elt1),AF.SecondaryStructure(elt2), "cityblock", "slow")]=(elt1,elt2)
    # keep=min(dict_dist.keys())
    # return dict_dist[keep],keep
    if len(l_struct2) >= len(l_struct1):
        res = pooling.starmap(calc_dist, [(elt2,l_struct1) for elt2 in l_struct2])
    else:
        res = pooling.starmap(calc_dist, [(elt1,l_struct2) for elt1 in l_struct1])
    
    pooling.terminate()
    
    res_fin=[]
    for elt in res:
        for elt1 in elt:
            res_fin.append(elt1)
    
    min=res_fin[0][0]
    for i,elt in enumerate(res_fin):
        if elt[0]<=min:
            choose=elt
            min=elt[0]
    return choose

def one_range_impact(struct_base, struct_test,aff_dash="min"):
    """Function designed to evaluate the impact of a single gap with the AptaMat distance
    
    struct_base and struct_test need to have a difference of 1 in their length.
    
    Returns an histogram of distances in regards to the position of the gap.
    
    aff_dash indicates where the dash is supposed to be inserted when displaying results.
    Can be: 'min'[default],'max','start','end'.
    """

    temp_struct="-"+struct_test
    L_struct_test=[temp_struct]
    for i in range(1,len(temp_struct)):
        temp_struct=del_str(temp_struct,i-1)
        temp_struct=insert_str(temp_struct,i,"-")
        L_struct_test.append(temp_struct)
        
    res=[]
    for elt in L_struct_test:
        res.append(AF.compute_distance_clustering(AF.SecondaryStructure(elt),AF.SecondaryStructure(struct_base), "cityblock", "slow"))

    aff_struct_test=[]
    if aff_dash!='end' and aff_dash!='start':
        i=0
        already=False
        while i != len(struct_test):
            if aff_dash=='min':
                if res[i] == min(res) and not already:
                    aff_struct_test.append("-")
                    already=True
                else:
                    aff_struct_test.append(struct_test[i])
                    i+=1
            elif aff_dash=='max':
                if res[i] == max(res) and not already:
                    aff_struct_test.append("-")
                    already=True
                else:
                    aff_struct_test.append(struct_test[i])
                    i+=1
    elif aff_dash == 'end':
        aff_struct_test=list(struct_test)
        aff_struct_test.append("-")
    elif aff_dash == 'start':
        aff_struct_test.append("-")
        for elt in struct_test:
            aff_struct_test.append(elt)
            
    fig, ax = plt.subplots()
    for i,elt in enumerate(list(struct_base)):
        ax.bar(i,res[i], color='blue')
        ax.text(i,-max(res)/20,elt)
        ax.text(i,-1.5*max(res)/20,aff_struct_test[i])
    ax.text(-len(struct_base)/10, -max(res)/20, "base")
    ax.text(-len(struct_base)/10, -1.5*max(res)/20, "compared")
    plt.show()
