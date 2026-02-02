import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(root_path)

import clustering_AptaMat as CL

structure_file="C:/Users/volni/Desktop/TX/AptaMat/clustering/datasets/triad/dataset_clustering_15x6_(cleaned_150x8reduced_nopk).dat"
cluster_dir = "C:/Users/volni/Desktop/TX/AptaMat/clustering/tests/tests_alignment/test_150x8/cleaned_nopseudoknots_reduced_6x15/clustering_results/distribution_results/"

structure_list,family=CL.initialize_dataset(structure_file)

def find_labels(cluster_file):
    with open(cluster_file) as f:
        lines = f.readlines()
    labels=[]
    families=[]
    ids=[]
    
    for i, line in enumerate(lines):
        if not (line.startswith("Optimal") or line.startswith("CLUSTER") or line=="\n"):
            l=line.split("   ")
            labels.append(l[0])
            families.append(l[1])
            ids.append(l[2])
        
    return labels, families, ids

for file in os.listdir(cluster_dir):
    if not file.endswith(".dat"):
        continue

    cluster_file = os.path.join(cluster_dir, file)
    print(f"Traitement de {file}")

    labels, families, ids = find_labels(cluster_file)
    dict_label = CL.build_label_dict(labels, family)
    
    CL.heatmap(family, labels, dict_label, title=file.replace("dataset_clustering_15x6_CLUSTERING_DISTRIBUTION_RESULTS_","Heatmap_color_").replace(".dat",".pdf"))
