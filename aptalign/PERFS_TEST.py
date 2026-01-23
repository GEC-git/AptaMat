import sys
import os

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','random_gen'))
sys.path.append(root_path)

import random_generator as RGEN
import AptAlign as AL
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np

def perf_test(max_length, var=10, depth=5):
    tab=[]
    for i in tqdm(range(10,max_length)):
        times=[]
        for tests in range(var):
            seq1,seq2=RGEN.double_struct_gen(i)
            t1=time.time()
            AL.clustering_opt_subdiv(seq1,seq2,depth)
            t2=time.time()
            tf=t2-t1
            times.append(tf)
        if len(times)==1:
            tab.append(times[0])
        else:       
            tab.append(times)
    return tab

lgt=1000
vr=10
X=[i for i in range(10,lgt)]
perfs=perf_test(lgt, var=vr)
Y=[np.mean(perfs[i]) for i in range(lgt-10)]
plt.plot(X,Y)

f_created=open("PERFRESULTS_var="+str(vr)+"_length="+str(lgt),'a')
f_created.write("RAW RESULTS\n")
f_created.write(str(perfs))
f_created.write("\nMEAN\n")
f_created.write(str(Y))
f_created.write("\nLENGTHS\n")
f_created.write(str(X))
f_created.close()

import csv
with open('perfs_test_d=5_u.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=",",
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for i,rows in enumerate(X):
        row=perfs[i]
        spamwriter.writerow(row)
        