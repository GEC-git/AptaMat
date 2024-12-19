import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptalign'))
sys.path.append(root_path)

import matplotlib.pyplot as plt
import AptAlign as AL

structure_file="/home/bcuvillier/Documents/AptaMat/clustering/tests/tests_align_nonaligned/Test4/CLEANED/data_clustering_BPRNA_RFAM_scraped_10000_CLUSTER_150x8_CLEANED_paired_aptaligned.dat"
structure_list, families = AL.initialize_dataset(structure_file)

fam_dict={}
for family in families:
    fam_dict[family]=[]
    
struct_obj_list=[]
for struct in structure_list:
    struct_obj_list.append(AL.Structure(struct[0],ident=str(struct[1]),fam=struct[2], AGU=struct[3]))

for struct in struct_obj_list:
    fam_dict[struct.family].append(struct)

plt.clf()
for fam in fam_dict.keys():
    dict_length={}
    for struct in fam_dict[fam]:
        if struct.length not in dict_length:
            dict_length[struct.length]=1
        else:
            dict_length[struct.length]+=1
    
    
    def getX(elt):
        return elt[0]
    
    xy=sorted(list(dict_length.items()), key=lambda pair : getX(pair))
    
    
    x=[]
    y=[]
    for elt in xy:
        x.append(elt[0])
        y.append(elt[1])
    plt.bar(x,y,width=1,label=fam.replace("_",""))



plt.title("Length distribution of all families")
plt.xlabel("Length")
plt.ylabel("Number of structures")
plt.legend()
plt.show()