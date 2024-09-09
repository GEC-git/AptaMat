# This is a chronological changelog for AptaMat

### WEEK 1 - 02/09/2024 -> 08/09/2024.

#### CHANGELOG

- Created AptaFast.py to separate the original algorithm from its optimised counter part.

- Added a compressed cache to store distances already calculated.
    
    The coordinates are hashed to an hexadecimal number then stored into a dictionary.

    Not very effective since the avoided calculation is just two substractions for manhattan distance and a square root for euclidean method.

- Added a small optimisation to ignore testing for common points in the two matrices, returning a distance equal to 0.

- Added a spiraling search inside the matrices to optimise searching for the nearest point.

    This method returns improper results when using big matrices. Small matrices are ok.

    ![Example](example_spiraling_search.png "Example of a search with datapoints extracted from the algorithm")

    Very effective in the cases of small matrices (~two to four time faster): very big matrices almost always take more time (~100 to 1000 times slower for 1000x1000)
    
    This is due to the number of points tested being dramatically larger the bigger the matrices are and eventually surpassing the number of points tested in older versions.
    
    From [wikipedia](https://en.wikipedia.org/wiki/Nearest_neighbor_search) : "The simplest solution to the NNS problem is to compute the distance from the query point to every other point in the database, keeping track of the "best so far". This algorithm, sometimes referred to as the naive approach, has a running time of O(dN), where N is the cardinality of S and d is the dimensionality of S. There are no search data structures to maintain, so the linear search has no space complexity beyond the storage of the database. Naive search can, on average, outperform space partitioning approaches on higher dimensional spaces."
        This is what's hapening here.
    
- Added this changelog.

#### FUTURE CHANGES AND IDEAS

- Optimising the search for the nearest point. 
    
    Determining if an approximate method is still satisfactory for big matrices.

    Using multiprocessing since the search for the nearest point for each point are completely independant.
    
    Determining the matrices of smallest distances for all points and using it directly in the calculations.
    
### WEEK 2 - 09/09/2024 -> 15/09/2024.

- Implemented multiprocessing to process all the nearest points calculations.
    - VERY effective in the case of bigger matrices. (720s to 40s on a 48 cores CPU with 1000x1000 matrices)
    
    - Now, the load is distributed across *n* cpu cores, *n* being chose by the user (with an educated prompt).
    
#### FUTURE CHANGES AND IDEAS

