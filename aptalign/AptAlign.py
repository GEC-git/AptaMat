import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptafast'))
sys.path.append(root_path)

import AptaFast as AF
#import time
import multiprocessing as mp
import matplotlib.pyplot as plt
#import numpy as np


### NO OOP METHOD

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
    """Function used to align two structures using a brute force method
    
    This calculates all possible alignments and returns the one with the smallest AptaMat distance.
    
    max_size is the size the structure will be aligned to.
        Is equal to the length of the biggest structure if under.
    """
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)
    print("Generating arrangement")
    l_struct1=arrangement(struct1,max_size)
    l_struct2=arrangement(struct2,max_size)
    print("S1: ",len(l_struct1),"| S2: ", len(l_struct2))
    print("This is a multiprocessed program, you have",mp.cpu_count(),"cores in your CPU.")
    nb=48#int(input("How much do you want to use? "))
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
    
    
def dynamic_one_range(struct, struct2, fill):  
    if fill != 0:
        print("Placing "+str(fill)+" gaps...")
        temp_struct=struct
        for j in range(fill):
            temp_struct="-"+temp_struct
            L_struct_test=[temp_struct]
            for i in range(1,len(temp_struct)):
                temp_struct=del_str(temp_struct,i-1)
                temp_struct=insert_str(temp_struct,i,"-")
                L_struct_test.append(temp_struct)
            res={}
            for i,elt in enumerate(L_struct_test):
                res[AF.compute_distance_clustering(AF.SecondaryStructure(elt),AF.SecondaryStructure(struct2), "cityblock", "slow")]=elt
                print("\r"+str(round((j/fill)*100,2))+"% | "+str(round(i/len(temp_struct)*100,2))+"%",end="\r")
            keep=min(res.keys())
            temp_struct=res[keep]
        struct = temp_struct
        return struct
    else:
        print("No gaps to place, continuing...")
        return struct
        
def dynamic_alignment(struct1,struct2,max_size=0):
    """
    Function to align two structures by adding gaps one by one in the optimal place.

    Parameters
    ----------
    struct1 : dotbracket
        First structure to align.
    struct2 : dotbracket
        Second structure to align.
    max_size : int, optional
        Size the structures take when aligned (controls the number of gaps). 
        The default is 0 then it is converted to match the biggest of the two structures.

    Returns
    -------
    The aligned structures with the resulting aptamat distance

    """
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)
    fill1=max_size-len(struct1)
    fill2=max_size-len(struct2)

    print("\nGenerating Alignment")
    print("\nFirst Structure")
    struct1_1=dynamic_one_range(struct1, struct2, fill1)
    
    print("\nSecond Structure")
    struct2_2=dynamic_one_range(struct2, struct1, fill2)
    
    print("\nFinished\n")
    
    return struct1_1, struct2_2, AF.compute_distance_clustering(AF.SecondaryStructure(struct1_1),AF.SecondaryStructure(struct2_2), "cityblock", "slow")



### NEW METHOD USING A STRUCTURAL ALPHABET AND OOP

def subdiv_finder(sequence,subdiv_param):
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
        
    return subdiv,new_seq

"""EXAMPLE:

`subdiv_finder(".((((((..((((.......)))).(((((.......))))).....(((.....)))))))))..",2)`

results in:
    subdiv = {6: '(',58: ')',5: '(',59: ')',4: '(',60: ')',3: '(',61: ')',2: '(',62: ')',1: '(',63: ')'},
    new_seq = '.######..((((.......)))).(((((.......))))).....(((.....)))######..'
"""

def slicer(sequence):
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
    for i in range(1,len(full_ranges)):
        dp_full_ranges.insert(i+insert,("SEP",[j for j in range(full_ranges[i-1][1][-1]+1,full_ranges[i][1][0])]))
        insert+=1
    if separator_end:
        dp_full_ranges.append(("SEP",[j for j in range(full_ranges[-1][1][-1]+1,len(sequence))]))
    
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
            
    
    
"""EXAMPLE:

`print(Structure(".((((((..((((.......)))).(((((.......))))).....(((.....))))))))).."))`

results in:
    
'''
Type: Structure
Length: 66
Raw Sequence: .((((((..((((.......)))).(((((.......))))).....(((.....)))))))))..
Subdiv sequence: .######..((((.......)))).(((((.......))))).....(((.....)))######..
Subdiv: {6: '(', 58: ')', 5: '(', 59: ')', 4: '(', 60: ')', 3: '(', 61: ')', 2: '(', 62: ')', 1: '(', 63: ')'}

Separators: 
Separator number 1:
Type: Separator
Length: 9
Order: 0
Sequence: .######..
Start: 0
Finish: 8
Subdiv index: 6

Separator number 2:
Type: Separator
Length: 1
Order: 2
Sequence: .
Start: 24
Finish: 24
Subdiv index: 0

Separator number 3:
Type: Separator
Length: 5
Order: 4
Sequence: .....
Start: 42
Finish: 46
Subdiv index: 0

Separator number 4:
Type: Separator
Length: 8
Order: 6
Sequence: ######..
Start: 58
Finish: 65
Subdiv index: 6

______________________________________

Patterns: 
Pattern number 1:
Type: Pattern
Length: 15
Order: 1
Start: 9
Finish: 23
Starting Sequence: ((((.......))))
Subdiv index: 0
Not yet aligned

Pattern number 2:
Type: Pattern
Length: 17
Order: 3
Start: 25
Finish: 41
Starting Sequence: (((((.......)))))
Subdiv index: 0
Not yet aligned

Pattern number 3:
Type: Pattern
Length: 11
Order: 5
Start: 47
Finish: 57
Starting Sequence: (((.....)))
Subdiv index: 0
Not yet aligned

______________________________________
Not yet aligned
'''

"""


class Structure():
    """
    Class representing an aligned or not dotbracket sequence.
    
    It is formed of two lists : one with its patterns and one with its separators.
    """
    def __init__(self, sequence):
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
        tbp+="Length: "+str(self.length)+"\n"
        tbp+="Raw Sequence: "+self.raw+"\n"
        tbp+="Subdiv sequence: "+str(self.raw_nosubdiv)+"\n"
        tbp+="Subdiv: "+str(self.subdiv)+"\n"
        tbp+="\nSeparators: \n"
        for i,elt in enumerate(self.separators):
            tbp+="Separator number "+str(i+1)+":\n"+str(elt)+"\n"
        tbp+="______________________________________\n"
        tbp+="\nPatterns: \n"
        for i,elt in enumerate(self.patterns):
            tbp+="Pattern number "+str(i+1)+":\n"+str(elt)+"\n"
        tbp+="______________________________________\n"

        if self.isaligned:
            tbp+="Aligned sequence: "+self.alignedsequence+"\n"
            tbp+="Aligned with:"+self.alignedwith.raw+"\n"
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
            
            
            for pattern in self.patterns():
                if not pattern.isaligned:
                    print("The patterns are not yet aligned, returning raw sequence")
                    return self.raw
            
            ordered=self.order_list()
            
            seqint=""
            for patsep in ordered:
                seqint+=patsep.sequence
            
            #assuming that subdiv will be updated after alignment.
            seq_dic={}
            for i,elt in enumerate(seqint):
                seq_dic[i]=elt
            for elt in self.subdiv.items():
                seq_dic[elt[0]]=elt[1]
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
        
        subdiv=0
        for elt in raw:
            if elt == "#":
                subdiv+=1
        
        self.subdiv_index=subdiv
        self.isaligned=False
        self.alignedsequence=raw
        
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
        
        subdiv=0
        for elt in raw:
            if elt == "#":
                subdiv+=1
        self.subdiv_index=subdiv
        
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
        self.alignedwith=acc_pat
        self.alignedsequence=new_seq
        self.length=len(new_seq)
        self.finish+=new_seq.count('-')
        self.isaligned=True
        
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
            tbp+="Aligned with: "+self.alignedwith.sequence+"\n"
        else:
            tbp+="Not yet aligned\n"
        return tbp

def add_gaps(current_order,added_gaps,order_list):
    
    for elt in order_list:
        if elt.nb > current_order:
            elt.start+=added_gaps
            elt.finish+=added_gaps

def pattern_alignment(struct1, struct2, pat1, pat2, order1, order2):
    """
    Used to align two patterns with dynamic alignment.
    """

    if len(pat1.sequence)==len(pat2.sequence):
        ms=len(pat1.raw)+1
    else:
        ms=0
        
    seq1,seq2,new_dist=dynamic_alignment(pat1.sequence, pat2.sequence,max_size=ms)
    
    pat1.aligned(pat2,seq1)
    pat2.aligned(pat1,seq2)
    
    added_gaps1=pat1.sequence.count("-")
    added_gaps2=pat2.sequence.count("-")
    
    if added_gaps1:
        add_gaps(pat1.nb, added_gaps1, order1)
        struct1.length+=added_gaps1
    if added_gaps2:
        add_gaps(pat2.nb, added_gaps2, order2)
        struct2.length+=added_gaps2
        
    print("Aligned:",pat1,pat2)

def matching_finder(struct1, struct2):
    """
    Used to pair up patterns when the number of patterns is not equal.
    """
    
def separator_compensating(struct1, struct2, matching):
    """
    Used after pattern aligning to add gaps in separators to match the final length and returning the aligned structures.
    
    """
    ordered1=struct1.order_list()
    ordered2=struct2.order_list()
    
    if struct1.length==struct2.length:
        
        for elt in matching.items():
            
            if not elt[0].start == elt[1].start:
                None
    
        
        
    #It is better to add as fewer gaps as possible.
    #We need to align patterns with each other in term of position in the structure.

def full_alignment(struct1, struct2):
    """
    
    Main function to calculate a pattern based alignment.
    
    struct1: Structure
    struct2: Structure
    
    ____________________
    METHOD:
        
        MATCHING
        
        - comparing number of patterns
        - comparing number of separators
        
        - comparing one by one the length and distance of the patterns and separators.
            - taking into account subdivisions with the subdiv_index.
        
        - Having a one by one match with all the patterns of the smallest structure.
            - A matched pattern in the bigger structure cannot have a smaller order than its match in the smaller structure.
            - If the number of pattern is the same in both structure, there is a single way to match each pattern.
            
        ALIGNING
        
        - Using dynamic alignment to align each pattern with its match.
            - mark them as aligned and input each other in the `self.alignedwith` variable.
        
        - When done, match the length of the two structures with gaps placed in the separators where it minimizes the distance.
        
        - Mark the structures as aligned and input each other in the `self.alignedwith` variable.
    ____________________
    
    """

    initial_dist=AF.compute_distance_clustering(AF.SecondaryStructure(struct1.sequence),AF.SecondaryStructure(struct2.sequence), "cityblock", "slow")
    
    pat_num1=len(struct1.patterns)
    pat_num2=len(struct2.patterns)
    
    # sep_num1=len(struct1.separators)
    # sep_num2=len(struct2.separators)
    
    matching={}
    
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
        
        print("Same structures, returning the same sequence")
        
        return struct1, struct2, initial_dist, initial_dist
        
    if pat_num1==pat_num2:
        print("Same number of patterns, aligning with the only match possible")
        for i,pat in enumerate(struct1.patterns):
            matching[pat]=struct2.patterns[i]
        
        for elt in matching.items():
            order1=struct1.order_list()
            order2=struct2.order_list()
            pattern_alignment(struct1,struct2,elt[0], elt[1],order1,order2)
            
        print("Compensating for length matching")
        separator_compensating(struct1, struct2)
        
        
        
        
