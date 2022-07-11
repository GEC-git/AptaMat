<img src="AptaMat.png" alt="AptaMat" width="200"/>
# AptaMAT

Purpose
-------------------

AptaMat is a simple script which aims to measure differences between DNA or RNA secondary structures. 
The method is based on the comparison of the matrices representing the two secondary structures to analyze, assimilable to dotplots. The dot-bracket
notation of the structure is converted in a half binary matrix showing width equal to structure's length.
Each matrix case (i,j) is filled with '1' if the nucleotide in position i is paired with the nucleotide in position j, with '0' otherwise. 

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

The `structures` argument must be a string formatted secondary structures. The first input structure is 
the template structure for the comparison. The following input are the compared structures. There are no input 
limitations. Quotes are necessary.

    usage: AptaMat.py structures [-h] "struct_1" "struct_2" ["struct_n" ...]

The `files` argument must be a formatted file. Multiple files can be parsed. The first structure encountered 
during the parsing is used as the template structure. The others are the compared structures.

    
    usage: AptaMat.py -files [-h] struct_file_1 [struct_file_n ...]


The input must be a text file, containing at least secondary structures, and accept additional 
information such as Title, Sequence or Structure index. If several files are provided, the function parses the files one
by one and always takes the first structure encountered as the template structure. Files must be formatted as follows: 


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
Compared structures need to have the same length as the Template structure.

For the moment, no features have been included to check whether the base pair is able to exist or not, according 
to literature. You must be careful about the sequence input and the base pairing associate.

The script accepts the extended dotbracket notation useful to compare pseudoknots or Tetrad. However, the resulting
distance might not be accurate.
