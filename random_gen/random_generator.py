from random import random


def random_gen(length = 1000,ratio_o_c=1):
    """
        Generates a random dotbracket notation of length *length* and density *ratio_o_c*
    
        What we call density is the number of brackets opened and closed.
        
        The density is dynamic so when a bracket is opened, a dot is more inclined to appear.
        When a dot is added, an open_bracket or a closing bracket is more likely to appear.
    """
    dot = "."
    open_bracket = "("
    closing_bracket= ")"

    opened=0
    output=""
    nb_open_bracket=0
    for i in range(length):
        if nb_open_bracket>=length-i:
            output+=closing_bracket
        else:
            r=random()
            if r > ratio_o_c:
                output+=dot
                ratio_o_c+=i/length**2
            elif r >= 1-ratio_o_c:
                output+=open_bracket
                nb_open_bracket+=1
                ratio_o_c-=(nb_open_bracket/length)
            else:
                if nb_open_bracket==0:
                    if random() >= 0.8:
                        output+=open_bracket
                        nb_open_bracket+=1
                        ratio_o_c-=(nb_open_bracket/length)
                    else:
                        output+=dot
                        ratio_o_c+=i/length**2
                else:
                    output+=closing_bracket
                    nb_open_bracket-=1
    return output

def dotbracket_verif(input_str):
    """
        This function verifies if all the opened brackets are actually closed.
        It appears that sometimes, the output of random_gen is invalid, we want to ignore those outputs.
    """
    nb=0
    for elt in input_str:
        if elt=="(":
            nb+=1
        if elt==")":
            nb-=1
            
    return not bool(nb)


def file_struct_gen(nb_struct=100):
    output=">Randomly generated structure file\n"
    """
>Test File Weight
--structure1--
.(((((..(.(((...))))...(((.....))).)))))
[ 0 ]
--structure2--
.(((((((..(((...)))..))(((.....))).)))))
[ 0.7 ]
--structure3--
.(((((......(....)......((.....))..)))))
[ 0.2 ]
--structure4--
(((.......)))......(((((((.....))).)))).
[ 0.1 ]
    """
    i=0
    while i < nb_struct:
        gen=random_gen()
        if dotbracket_verif(gen):
            output+="--structure"+str(i+1)+"--\n"
            output+=gen+"\n"
            output+="[ "+str(random())+" ]\n"
            i+=1
    return output

f=open("RandomWeightFile.fa","a")
f.write(file_struct_gen())
f.close()
print("Successfuly created a randomly generated weighted structure file.")
        