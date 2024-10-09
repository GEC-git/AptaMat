import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptafast'))
sys.path.append(root_path)

import numpy as np
import AptaFast as AF

import time
import multiprocessing


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


def double_mode_brute_force(struct1,struct2,max_size):
    
    if len(struct1)>max_size: max_size=len(struct1)
    if len(struct2)>max_size: max_size=len(struct2)
    
    fill1=max_size-len(struct1)
    fill2=max_size-len(struct2)
    
    test="--.((..))...(((..))).."
    t=   "---------.(..(((())))."
    
    