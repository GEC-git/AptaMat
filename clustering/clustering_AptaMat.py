import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptamat2.0'))
sys.path.append(root_path)

import numpy as np
import pandas as pd
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import calinski_harabasz_score, silhouette_score, adjusted_rand_score
import AptaMat2 as AF
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import time
import multiprocessing
from vispy import scene
from vispy import app
import argparse




def build_label_dict(labels, family):
    """
    Use unorganised labels list to build a dictionary of family asssociated labels

    Parameters
    ----------
    labels :  list
        Input list of labels
    family : list
        Input list of family id

    Returns
    -------
        label_dict : dict
            Updated list of labels
    """
    family_range = []
    family_name = []
    sum_value = 0
    for n, (key, value) in enumerate(family.items()):
        family_name.append(key)
        family_range.append(sum_value)
        sum_value += value

    label_dict = {}
    for n, range_fami in enumerate(family_range):
        newdict = {}
        try:
            set(labels[family_range[n]:family_range[n + 1]])
        except:
            for i in set(labels[family_range[n]:]):
                newdict[i] = labels[family_range[n]:].count(i)
            label_dict[family_name[n]] = newdict
        else:
            for i in set(labels[family_range[n]:family_range[n + 1]]):
                newdict[i] = labels[family_range[n]:family_range[n + 1]].count(i)
            label_dict[family_name[n]] = newdict

    return label_dict


def renumber_by_rank(labels):
    """
    Renumber labels in dictionary in prevailing order

    Parameters
    ----------
    labels : list
        uncorrected label list

    Returns
    -------
        new_labels : list
            corrected label list
    """
    dcc = {}
    for l in labels:
        dcc[l] = dcc.get(l, 0) + 1

    ranked = sorted(dcc.items(), key=lambda x: x[1], reverse=True)

    new_labels = []
    for element in labels:
        for i, r in enumerate(ranked):
            if element == r[0]:
                new_labels.append(i)     
    return new_labels


def affinity_calculation(distance_matrix, standard_labels, sigma, depth):
    """Used for parallelizing the affinity_propagation calculation"""
    affinity_matrix = np.exp(- distance_matrix ** 2 / (2. * sigma ** 2))
    clustering = AffinityPropagation(damping=0.52, affinity="precomputed", convergence_iter=10, max_iter=depth,
                                     random_state=0).fit(affinity_matrix)
    
    #sub_aff_prop.append(clustering)
    calinski = calinski_harabasz_score(distance_matrix, clustering.labels_)
    silhouette = silhouette_score(affinity_matrix, clustering.labels_)
    acc_score = adjusted_rand_score(standard_labels, clustering.labels_)
    
    return (clustering, calinski, silhouette, acc_score, sigma, affinity_matrix)

def optimised_affinity_propagation(distance_matrix, CORE, depth, standard=None, sigma=np.arange(1, 10, 0.1)):
    ### Acquire expected clustering labels
    if standard is not None:
        standard_labels = []
        for v in standard:
            standard_labels.append(v)
    else:
        standard_labels=[i for i in range(len(distance_matrix))]

    if isinstance(sigma, (list, tuple, np.ndarray)):
        sigma_iter = sigma
    elif isinstance(sigma, (int, float)):
        sigma_iter = [sigma]

    ### AffinityPropagation cluster creation applying various sigma value
    ### Each iteration results in calculation of various clustering quality metrics
    sub_aff_prop = []
    aff_prop_calinski_best = 0
    silhouette_best = -1
    aff_prop_clust_best = None
    acc_best = 0
    sigma_best = 0
    print("Recreating pool for affinity calculation.")
    pool = multiprocessing.Pool(CORE)
    results=[]
    print("Job started")
    for result in pool.starmap(affinity_calculation,
                               [(distance_matrix, standard_labels, sigma, depth) for sigma in sigma_iter]):
        results.append(result)
    pool.terminate()
    print("Job Finished, fetching the best result.")
    for elt in results:
        
        acc_score=elt[3]
        calinski=elt[1]
        silhouette=elt[2]
        if acc_best < acc_score:
            acc_best = acc_score
            aff_prop_calinski_best = calinski
            silhouette_best = silhouette
            aff_prop_clust_best = elt[0]
            sigma_best = elt[4]
            aff=elt[5]

        if aff_prop_calinski_best < calinski and silhouette_best < silhouette:
            aff_prop_calinski_best = calinski
            silhouette_best = silhouette
            aff_prop_clust_best = elt[0]
            sigma_best = elt[4]
            aff=elt[5]

        elif aff_prop_calinski_best > calinski and silhouette_best < silhouette:
            if calinski + (5 * aff_prop_calinski_best / 100) > aff_prop_calinski_best:
                aff_prop_calinski_best = calinski
                silhouette_best = silhouette
                aff_prop_clust_best = elt[0]
                sigma_best = elt[4]
                aff=elt[5]

        elif aff_prop_calinski_best < calinski and silhouette_best > silhouette:
            if silhouette + (20 * silhouette_best / 100) > silhouette_best:
                aff_prop_calinski_best = calinski
                silhouette_best = silhouette
                aff_prop_clust_best = elt[0]
                sigma_best = elt[4]
                aff=elt[5]
                
    return aff, aff_prop_clust_best, aff_prop_calinski_best, silhouette_best, acc_best, sigma_best, sub_aff_prop


### Initialize dataset
def initialize_dataset(structure_file):
    structure_list = []
    family = {}
    with open(structure_file, 'r') as file:
        for line in file:
            content = line.strip().split()
            #print(content)
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
                                                           id=content[1].split('.')[0],fam=str(content[0]))
                    # AptaMat._create_fasta(structure)
                    structure_list.append(structure)
    return structure_list,family

def calculation(structure_list, CORE, speed, depth, sigma_range):
    
    ### N for matrix size
    N = len(structure_list)
    
    ### Calculate AptaMat distance for each
    ### Each structure comparison result in a tuple (struct1, struct2, AptaMat distance)
    
    #AF.compute_distance(struct_1, struct_2, method, nb_pool, pool, speed)
    start = time.time()
    print("Job started",time.asctime())
    results = []
    pool = multiprocessing.Pool(CORE)
    for result in pool.starmap(AF.compute_distance_clustering,
                               [(struct1, struct2,"cityblock",speed) for struct1 in structure_list for struct2 in structure_list]):
        results.append(result)
    
    pool.terminate()
    
    ### Build distance matrix using AptaMat distance in 'results' tuples
    matrix_element = []
    for i in results:
        if i == None:
            matrix_element.append(0)
        else:
            matrix_element.append(float(i))
        
    dist_matrix = np.array(matrix_element).reshape(N, N)
    
    
    ### Acquire data from Affinity Propagation clustering
    affinity_matrix, aff_prop_clust_best, aff_prop_calinski_best, silhouette_best, acc_best, sigma_best, sub_aff_prop = \
        optimised_affinity_propagation(dist_matrix, CORE, depth, sigma=np.arange(1, sigma_range, 0.1)) 
     
    end = time.time()
    print("Job finished",time.asctime())
    print("Time elapsed = ", time.strftime("%H:%M:%S", time.gmtime(end-start)))
    return affinity_matrix, aff_prop_clust_best, aff_prop_calinski_best, silhouette_best, acc_best, sigma_best, sub_aff_prop


def affinity_visualization_GPU(affinity_matrix, structure_list):
    precision=len(structure_list)
    tristogram = np.zeros((precision, precision, precision), dtype=np.uint8)
    for i,elt in enumerate(affinity_matrix):
        for j,num in enumerate(elt):
            if i<=j:
                tristogram[int(num*(precision-1)),i,j]=num*100
            else:
                tristogram[int(num*(precision-1)),i,j]=0
                
    canvas = scene.SceneCanvas(keys='interactive', bgcolor='w')
    
    view = canvas.central_widget.add_view()
    volume = scene.visuals.Volume(tristogram, parent=view.scene,cmap="RdBu")

    view.camera = scene.cameras.TurntableCamera(parent=view.scene,up='z', fov=60)
    
    yax = scene.Axis(pos=[[-25, 0], [-25, precision]], tick_direction=(-1, 0),
                     font_size=precision*10, axis_color='k', tick_color='k', text_color='k',
                     parent=view.scene, axis_label="Y",axis_label_margin=precision*20,
                     major_tick_length=precision*5,minor_tick_length=precision*2,tick_label_margin=precision*10)
    
    xax = scene.Axis(pos=[[0, -25], [precision, -25]], tick_direction=(0,-1),
                     font_size=precision*10, axis_color='k', tick_color='k', text_color='k',
                     parent=view.scene, axis_label="X",axis_label_margin=precision*20,
                     major_tick_length=precision*5,minor_tick_length=precision*2,tick_label_margin=precision*10)
    
    if __name__ == '__main__':
        canvas.show()
        app.run()


def affinity_visualization_CPU(affinity_matrix,structure_list):
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x=np.arange(0,len(structure_list),1)
    y=np.arange(0,len(structure_list),1)
    hist, xedges, yedges = np.histogram2d(x, y, bins=4, range=[[0, len(structure_list)], [0, len(structure_list)]])
    zdir=(0,-1,0)
    ydir=(1,0,0)
    for i,elt in enumerate(structure_list):
        ax.text(i, -50, 0, elt.id, zdir,size=3)
        ax.text(len(structure_list)+50,i,0,elt.id,ydir,size=3)
    
    
    xpos, ypos = np.meshgrid(x,y, indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    
    # Construct arrays with the dimensions for the 16 bars.
    dx = dy = 0.5 * np.ones_like(zpos)
    dz = np.array([affinity_matrix[i][j] if i>=j else 0 for i in range(len(structure_list)) for j in range(len(structure_list))])
    offset = dz + np.abs(dz.min())
    fracs = offset.astype(float)/offset.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.jet(norm(fracs.tolist()))
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=color_values)
    
    plt.show()

def heatmap(family, labels, dict_label):
### Heatmap setup
    #df = pd.DataFrame(family, index=list(set(labels)), columns=dict_label.keys())
    #print(dict_label)
    df=pd.DataFrame(dict_label)
    s = df.sum()
    df = df[s.sort_values(ascending=False).index[:]]
    family_np = df.to_numpy()
    family_np=np.nan_to_num(family_np)
    
    family_percent = family_np / family_np.sum(axis=0) * 100
    family_np_t = np.transpose(family_np)
    family_p_t = np.transpose(family_percent)
    clean_labels = [i.replace('_', ' ') for i in df.columns]
    
    binary_m = plt.get_cmap('jet')
    colormap = ListedColormap(binary_m(np.linspace(0.2, 1, 100)))
    colormap.set_under(color='white')
    fig, ax = plt.subplots()
    im = ax.imshow(family_p_t, cmap=colormap, vmin=0.9)
    ax.set_yticks(np.arange(len(df.columns)), labels=clean_labels)
    # ax.set_yticks(np.arange(len(rfam)), labels=rfam)
    ax.set_xticks(np.arange(0, len(df)), labels=list(np.arange(0, len(df))))
    #print(df)
    for i in range(len(df.columns)):
        for j in np.arange(0, len(df)):
            if round(family_p_t[i, j]) == 0:
                pass
            elif round(family_p_t[i, j]) < 75:
                text = ax.text(j, i, int(family_np_t[i, j]),
                               ha="center", va="center")
            else:
                text = ax.text(j, i, int(family_np_t[i, j]),
                               ha="center", va="center", color="w")
    
    plt.text(len(df)+1, 1, 'Occupancy (%)', fontsize=10)
    cax = plt.axes([0.98, 0.295, 0.02, 0.4])
    plt.colorbar(im, cax)
    fig.savefig('HeatMap_Color.pdf', dpi=600, bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser(description="This clustering algorithm uses AptaFast to determine"
                                     "a distribution of structures inputed from a file.")
    parser.add_argument('-fp',
                        '--filepath',
                        type=str,
                        nargs='+',
                        help='Input file containing structures to be clustered.')
    
    parser.add_argument('-speed', 
                        default='slow',
                        help="Using greedy or non greedy depth calculation",
                        nargs='?',
                        choices=['slow','quick'])
    
    parser.add_argument('-visu',
                        default=None,
                        help="Using GPU or CPU accelerated affinity_matrix visualization. Default : None",
                        nargs='?',
                        choices=['GPU','CPU'])
    
    parser.add_argument('-cv',
                        '--cluster_visualization',
                        help="Displaying heatmap at the end",
                        action="store_true")
    
    parser.add_argument('-d',
                        '--depth',
                        type=int,
                        default=1000,
                        nargs='+',
                        help="Depth of clustering calculation.")
    
    parser.add_argument('-sr',
                        '--sigma_range',
                        type=int,
                        default=10,
                        nargs='+',
                        help="Range of clustering calculation.")
    
    args = parser.parse_args()
    if isinstance(args.depth,list):
        depth=args.depth[0]
    else:
        depth=args.depth
        
    if isinstance(args.sigma_range,list):
        sigma_range=args.sigma_range[0]
    else:
        sigma_range=args.sigma_range
    ### Structure file to be used
    structure_file=""
    for elt in args.filepath:
        structure_file+=str(elt)
    
    print("This is a multiprocessed algorithm.")
    print("You have",multiprocessing.cpu_count(),"cores in your CPU.")
    CORE = int(input("Please input the number of cores you want to use:"))
    
    structure_list,family=initialize_dataset(structure_file)
    
    affinity_matrix, aff_prop_clust_best, aff_prop_calinski_best, silhouette_best, acc_best, sigma_best, sub_aff_prop=calculation(structure_list, CORE, args.speed, depth, sigma_range)
    
    
    ### Print Optimal values obtained from affinity propagation clustering
    print('Optimal Calinski Harabasz index =', aff_prop_calinski_best)
    print("Optimal Silhouette score =", silhouette_best)
    print('Optimal Sigma =', sigma_best)
    labels = aff_prop_clust_best.labels_
    tbw="CLUSTER   FAMILY   ID   DOTBRACKET\n"
    for i,struct in enumerate(structure_list):
        tbw+=str(labels[i])+"   "+struct.family+"   "+struct.id+"   "+struct.dotbracket+"\n"
    
    f_created=open(structure_file.replace(".dat","")+"_CLUSTERING_DISTRIBUTION_RESULTS.dat",'a')
    f_created.write(tbw)
    f_created.close()
    labels = renumber_by_rank(labels)
    dict_label = build_label_dict(list(labels),family)
    
    if args.visu is not None:
        if args.visu=="GPU":
            affinity_visualization_GPU(affinity_matrix, structure_list)
        if args.visu=="CPU":
            affinity_visualization_CPU(affinity_matrix, structure_list)
    
    if args.cluster_visualization:
        heatmap(family, labels, dict_label)
    
if __name__ == '__main__':
    main()
