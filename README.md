# AptaMAT

Purpose
-------------------

AptaMat is a simple script aimed to measure differences between DNA or RNA secondary structures. 
The method is inspired by dotplot approach developed in Ivry's paper (Ivry et al., 2009). The dotbracket
notation of the structure is converted in a half binary matrix showing width equal to structure's length.
From "0" and "1" we are able to recover base paired position in structures. 

The differences between matrices is calculated by applying Manhattan distance on each point in the template matrix 
against all the points from the compared matrix. This calculation is repeated between compared matrix and template
matrix to handle all the differences. Both calculation are then sum up and divided by the sum of all the points in 
both matrices.

Dependencies
------------

AptaMat have been written in Python 3.8+

Two Python modules are needed :

- [NumPy](https://numpy.org/)
- [scipy](https://www.scipy.org/)

These can be installed by typing in the command prompt either :

    ./setup
or

    pip install numpy
    pip install scipy

Use of [Anaconda](https://docs.conda.io/en/latest/#) is highly recommended.

Usage
------------

AptaMat is a flexible Python script which can take several arguments:

- `structures` followed by secondary structures written in dotbracket format
- `files` followed by path to formatted files containing one, or several secondary structures in dotbracket format

Both `structures` and `files` are independent functions in the script and cannot be called at the same time.

    usage: AptaMAT.py [-h] [-structures STRUCTURES [STRUCTURES ...]] [-files FILES [FILES ...]] 

The `structures` argument must take as argument string formatted secondary structures. The first input structure is 
the Template structure for the comparison. The following input are the Compared structures. There are no input 
limitation. Quotes are necessary.

    usage: AptaMat.py structures [-h] "struct_1" "struct_2" ["struct_n" ...]

The `files` argument take as argument formatted file. Multiple file can be parsed. The first structure encountered 
during the parsing is the Template structure. The others are the Compared structures.

    
    usage: AptaMat.py -files [-h] struct_file_1 [struct_file_n ...]


The input must be a text file, containing at least secondary structures, and accept additional 
information such as Title, Sequence or Structure index. If several files are provided, the function parse the files one
by one and always take the first structure encountered as the Template structure. Files must be formatted as follows: 


    >5HRU
    TCGATTGGATTGTGCCGGAAGTGCTGGCTCGA
    --Template--
    ((((.........(((((.....)))))))))
    --Compared--
    .........(((.(((((.....))))).)))

Examples
------------

### structures function
First introducing a simple example with 2 structures:

    $ AptaMat.py -structures "(((...)))" "((.....))"
     (((...)))
     ((.....))
    > AptaMat : 0.08
    
Then, it is possible to input several structures:
    
    $ AptaMat.py -structures "(((...)))" "((.....))" ".(.....)." "(.......)"
     (((...)))
     ((.....))
    > AptaMat : 0.08
    
     (((...)))
     .(.....).
    > AptaMat : 0.2
    
     (((...)))
     (.......)
    > AptaMat : 0.3

### files function
Taking the above file example:

    $ AptaMat.py -files example.fa
    5HRU
    Template - Compared
     ((((.........(((((.....)))))))))
     .........(((.(((((.....))))).)))
    > AptaMat : 0.1134453781512605

Note
------------
Compared structures have to show the same length as the Template structure.

For the moment, no features have been included to check whether the base pair is able to exist or not, according 
to literature. You must be careful about the sequence input and the base pairing associate.

The script accepts the extended dotbracket notation useful to compare pseudoknots or Tetrad. However, the resulting
distance might not be accurate.
