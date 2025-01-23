import subprocess
import os


def locarna_pairwise(struct1,struct2):
    tbw1=">"+struct1.id+"\n"
    tbw1+=struct1.sequence+"\n"
    tbw1+=struct1.dotbracket
    f_created=open("temp_struct1"+struct1.id+str(os.getpid())+".fa",'a')
    f_created.write(tbw1)
    f_created.close()
    tbw2=">"+struct2.id+"\n"
    tbw2+=struct2.sequence+"\n"
    tbw2+=struct2.dotbracket
    f_created=open("temp_struct2"+struct2.id+str(os.getpid())+".fa",'a')
    f_created.write(tbw1)
    f_created.close()
    
    ident1, sequence1al, dotbracket1al, ident2, sequence2al, dotbracket2al = locarna_from_file("temp_struct1"+struct1.id+str(os.getpid())+".fa", "temp_struct2"+struct2.id+str(os.getpid())+".fa")
    
    os.remove("temp_struct1"+struct1.id+str(os.getpid())+".fa")
    os.remove("temp_struct2"+struct2.id+str(os.getpid())+".fa")
    
    return dotbracket1al, dotbracket2al

def locarna_from_file(file1, file2):
    result = subprocess.run(["locarna",file1,file2], capture_output=True, text=True)
    lines=result.stdout.splitlines()
    del lines[0]
    del lines[0]
    del lines[2]
    del lines[-1]
    num=1
    for line in lines:
        if num == 1:
            j=0
            char=line[j]
            ident1=""
            while char != " ":
                ident1+=char
                char=line[j]
                j+=1
            while char == " ":
                char=line[j]
                j+=1
            sequence1=line[j:]
            num+=1
        elif num==2:
            j=0
            char=line[j]
            ident2=""
            while char != " ":
                ident2+=char
                char=line[j]
                j+=1
            while char == " ":
                char=line[j]
                j+=1
            sequence2=line[j:]
            num+=1
        elif num==3:
            j=0
            char=line[j]
            ident1=""
            while char != " ":
                ident1+=char
                char=line[j]
                j+=1
            while char == " ":
                char=line[j]
                j+=1
            dotbracket1=line[j:]
            num+=1
        elif num==4:
           j=0
           char=line[j]
           ident2=""
           while char != " ":
               ident2+=char
               char=line[j]
               j+=1
           while char == " ":
               char=line[j]
               j+=1
           dotbracket2=line[j:]
           num+=1
    
    return ident1, sequence1, dotbracket1, ident2, sequence2, dotbracket2
