import subprocess
import os
import sys

current_dir = os.path.dirname(__file__)
root_path = os.path.abspath(os.path.join(current_dir, '..','aptamat2.0'))
sys.path.append(root_path)

import AptaMat2 as AF

def locarna_pairwise(struct1,struct2):

    tbw1=">"+str(struct1.id)+"\n"
    tbw1+=struct1.dotbracket
    f_created=open("temp_struct1"+str(struct1.id)+str(os.getpid())+".fa",'a')
    f_created.write(tbw1)
    f_created.close()
    tbw2=">"+str(struct2.id)+"\n"
    tbw2+=struct2.dotbracket
    f_created=open("temp_struct2"+str(struct2.id)+str(os.getpid())+".fa",'a')
    f_created.write(tbw2)
    f_created.close()
    
    ident1, dotbracket1al, ident2, dotbracket2al = locarna_from_file("temp_struct1"+str(struct1.id)+str(os.getpid())+".fa", "temp_struct2"+str(struct2.id)+str(os.getpid())+".fa")
    
    os.remove("temp_struct1"+str(struct1.id)+str(os.getpid())+".fa")
    os.remove("temp_struct2"+str(struct2.id)+str(os.getpid())+".fa")
    
    return dotbracket1al, dotbracket2al

def locarna_from_file(file1, file2):
    result=subprocess.run(["locarna","--pp="+file1.replace(".fa","")+file2.replace(".fa","")+".pp",file1,file2],capture_output=True,text=True)
    del result
    lines=open(file1.replace(".fa","")+file2.replace(".fa","")+".pp").readlines()
    ### GAPS ARE ALSO IN SEQUENCE SO TO REVAMP
    for i,line in enumerate(lines):
        line=line.replace("\n","")
        if not(line.startswith("#") and len(line) == 0):
            if i==2:
                ident1=""
                j=2
                while line[j]!=" ":
                    ident1+=line[j]
                    j+=1
                while line[j]==" ":
                    j+=1
                dotbracket1=line[j:]
            if i == 3:
                ident2=""
                j=2
                while line[j]!=" ":
                    ident2+=line[j]
                    j+=1
                while line[j]==" ":
                    j+=1
                dotbracket2=line[j:]
    
    os.remove(file1.replace(".fa","")+file2.replace(".fa","")+".pp")
    return ident1, dotbracket1, ident2, dotbracket2
