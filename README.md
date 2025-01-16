<img src="triptych.png" alt="AptaMat2.0" height="300"/>

# Purpose
-------------------

This repository host the AptaMat2.0, AptAlign and clustering algorithms based on the first version of AptaMat.
You can also find other scripts to convert datasets or generating ones randomly.

____

**AptaMat2.0** is a simple optimised script which aims to measure differences between DNA or RNA secondary structures. 
The method is based on the comparison of the matrices representing the two secondary structures to analyze, assimilable to dotplots. 
The dot-bracket notation of the structure is converted in a half binary matrix showing width equal to structure's length.
Each matrix case (i,j) is filled with '1' if the nucleotide in position i is paired with the nucleotide in position j, with '0' otherwise. 

The differences between matrices is calculated by applying Manhattan distance on each point in the template matrix 
against all the points from the compared matrix. This calculation is repeated between compared matrix and template
matrix to handle all the differences. Both calculation are then sum up with the number of gaps encountered and divided 
by the sum of all the points in both matrices.

AptaMat can handle extended dot-bracket notation and every additional bracket is converted into coordinates for the matrix.

AptaMat can also compare structures of different length. However, we recommend to work with structure of same length. Our 
The algorithm includes gap understanding, where each gap is considered as an additional penalized unpaired nucleotide.

____

**AptAlign** is an alignement algorithm based on pattern recognition inside of the dotbracket notation.
It is designed to minimize the AptaMat distance between two DNA or RNA secondary structures by inserting gaps.

This algorithm slices the dotbracket notation into notable patterns for each structure, then match the patterns between the two with three passes of calculation using a scoring matrix.
It then aligns the patterns inside of the structures and return aligned dotbracket notations for each structures.

AptAlign can align two structures or an ensemble of structures inside of a CLUSTER file.

____


The **clustering** algorithm can be used with AptAlign to generate clusters from a dataset.
It uses the AptaMat distance as a similarity index for the affinity propagation clustering provided by *scikit-learn*.

It can also display the affinity matrix when finished using a CPU or GPU accelerated 3D matrix visualizer.

This algorithm then gives results inside of a CLUSTER file or as an heatmap if a quality test is performed using a dataset with known families.

# Dependencies
------------

AptaMat2.0, AptAlign and clustering have been written in Python 3.8+

Those Python modules are needed for all algorithms:

- [NumPy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [Multiprocessing](https://docs.python.org/3/library/multiprocessing.html)
- [varnaapi](https://amibio.gitlabpages.inria.fr/varna-api/)
- [pandas](https://pandas.pydata.org/)
- [scikit-learn](https://scikit-learn.org/stable/)

Add for clustering only:
- [matplotlib](https://matplotlib.org/)
- [vispy](https://vispy.org/)

Use of [Anaconda](https://docs.conda.io/en/latest/#) is highly recommended.

# Usage
------------

## AptaMat2.0

AptaMat2.0 is a flexible Python script which can take several arguments:

- `-structures` followed by secondary structures written in dotbracket format
- `-weigths` (Optionnal) followed by weight values between 0 and 1 to indicate optionnal weight indices
- `-files` followed by path to formatted files containing one, or several secondary structures in dotbracket format
- `-ensemble` (Optionnal) which indicates whether the input secondary structures are part of an ensemble
- `-method` indicates the spatial distance method choose for AptaMat, by default cityblock and alternatively euclidean
- `-speed` indicates the risk taken by the algorithm when calculating the searchg depth. (default: slow) Can be set to quick if the user is confident in its data.

      usage: AptaMat2.py [-h] [-v] [-speed [{slow,quick}]] [-structures STRUCTURES [STRUCTURES ...]] [-weights WEIGHTS [WEIGHTS ...]] [-files FILES [FILES ...]] [-ensemble] [-method [{cityblock,euclidean}]]
      
Both `structures` and `files` are independent functions in the script and cannot be called at the same time.

The `structures` argument must be a string formatted secondary structures array. The first input structure is 
the template structure for the comparison. The following input are the compared structures. There are no input 
limitations. Quotes are necessary.


      usage: AptaMat2.py -structures STRUCTURES [STRUCTURES ...]


The `weight` optionnal argument must be an array of float in 0 to 1 range showing identical size than input `structures` array. 
This argument is not compatible with `files` as the script is expecting this information to be in the input file. 


      usage: AptaMat2.py -structures STRUCTURES [STRUCTURES ...] -weigths WEIGHTS [WEIGHTS ...]
    
    
The `files` argument must be a formatted file. Multiple files can be parsed. The first structure encountered 
during the parsing is used as the template structure. The others are the compared structures.

    
      usage: AptaMat.py -files FILES [FILES ...]
    

The input must be a text file, containing at least secondary structures, and accept additional 
information such as Title, Sequence, Structure index and Weight. If several files are provided, the function parses the files one
by one and always takes the first structure encountered as the template structure. Files must be formatted as follows: 

      >5HRU
      TCGATTGGATTGTGCCGGAAGTGCTGGCTCGA
      --Template--
      ((((.........(((((.....)))))))))
      [ weight ]
      --Compared--
      .........(((.(((((.....))))).)))
      [ weight ]
      ..........((.((((.......)))).)).
      [ weight ]


`ensemble` is an optionnal argument which allow to calculate AptaMat distance value for an ensemble of structure
instead of calculating pairwise distance.


      usage: AptaMat.py -structures STRUCTURES [STRUCTURES ...] -weigths WEIGHTS [WEIGHTS ...] -ensemble
          or
      usage: AptaMat.py -files FILES [FILES ...] -ensemble

# Note
------------

For the moment, no features have been included to check whether the base pair is able to exist or not, according 
to literature. You must be careful about the sequence input and the base pairing associated.


# Citation
------------
If you are using AptaMat in your research, please support us by citing us : Thomas Binet, Bérangère Avalle, Miraine Dávila Felipe, Irene Maffucci, AptaMat: a matrix-based algorithm to compare single-stranded oligonucleotides secondary structures, Bioinformatics, Volume 39, Issue 1, January 2023, btac752, https://doi.org/10.1093/bioinformatics/btac752

If you are using AptaMat2.0, an article is on the way.