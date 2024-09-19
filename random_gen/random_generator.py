from random import random


def random_gen(length=1000, ratio_o_c=0.8, bias=0.6):
    """
        Generates a random dotbracket notation of length *length*, density *ratio_o_c* and bias *bias*.
    
        What we call density is the number of brackets opened and closed.
        
        The density is dynamic so when a bracket is opened, a dot is more inclined to appear.
        When a dot is added, an open_bracket or a closing bracket is more likely to appear.
        
        The case `()` is not possible to appear for realistic purposes.
    """
    dot = "."
    open_bracket = "("
    closing_bracket= ")"
    just_opened=False
    output=""
    nb_open_bracket=0
    i=0
    while i != length:
        if length-i==2 or length-i==1 or length-i==0:
            if nb_open_bracket==0:
                output+=dot
                ratio_o_c+=i/length**2
                just_opened=False
                i+=1
            else:
                if just_opened:
                    output+=dot
                    ratio_o_c+=i/length**2
                    just_opened=False
                    i+=1
                else:
                    if nb_open_bracket>=length-i:
                        output+=closing_bracket
                        nb_open_bracket-=1
                        just_opened=False
                        i+=1
                    else:
                        if random() < bias:
                            output+=dot
                            ratio_o_c+=i/length**2
                            just_opened=False
                            i+=1
                        else:
                            output+=closing_bracket
                            nb_open_bracket-=1
                            just_opened=False
                            i+=1
        else:
            if nb_open_bracket>=length-i and not(just_opened):
                output+=closing_bracket
                nb_open_bracket-=1
                just_opened=False
                i+=1
            elif nb_open_bracket>=length-i and just_opened:
                output+=dot
                ratio_o_c+=i/length**2
                just_opened=False
                i+=1
            else:
                r=random()
                if r > ratio_o_c:
                    output+=dot
                    ratio_o_c+=i/length**2
                    just_opened=False
                    i+=1
                elif r >= 1-ratio_o_c:
                    output+=open_bracket
                    nb_open_bracket+=1
                    just_opened=True
                    ratio_o_c-=(nb_open_bracket/length)
                    i+=1
                else:
                    if nb_open_bracket==0:
                        if random() >= bias:
                            output+=open_bracket
                            nb_open_bracket+=1
                            ratio_o_c-=(nb_open_bracket/length)
                            just_opened=True
                            i+=1
                        else:
                            output+=dot
                            ratio_o_c+=i/length**2
                            just_opened=False
                            i+=1
                    elif not(just_opened):
                        output+=closing_bracket
                        nb_open_bracket-=1
                        just_opened=False
                        i+=1
                    
    return output

def dotbracket_verif(input_str,length):
    """
        This function verifies if all the opened brackets are actually closed.
        It appears that sometimes, the output of random_gen is invalid, we want to ignore those outputs.
    """
    nb=0
    if len(input_str)!=length:
        return False
    for elt in input_str:
        if elt=="(":
            nb+=1
        if elt==")":
            nb-=1
            
    return not bool(nb)


def file_struct_gen(nb_struct,lgth):
    output=">Randomly generated structure file\n"
    i=0
    while i < nb_struct:
        gen=random_gen(length=lgth)
        if dotbracket_verif(gen,lgth):
            output+="--structure"+str(i+1)+"--\n"
            output+=gen+"\n"
            output+="[ "+str(random())+" ]\n"
            i+=1
            
    return output
lgth=1000
nb_struct=10000
f=open("RandomWeightFile"+str(nb_struct)+"-"+str(lgth)+".fa","a")
f.write(file_struct_gen(nb_struct,lgth))
f.close()
print("Successfuly created a randomly generated weighted structure file.")
        