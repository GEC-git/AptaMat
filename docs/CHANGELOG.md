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

- Implemented double list search in place of spiraling search.

    - This type of search uses the properties of naive search on sub matrices centered around the originating point. Since the points are very dispersed in the matrices, there's a good chance to find the closest one near the originating point.
    
    - With multiprocessing, it is around 2x faster for small matrices (~20x20 and a pool of 2 cores) and also 2x faster for bigger matrices (~1000x1000 and a pool of 4 cores).
    
    - For bigger matrices, this method still returns wrong results but for smaller ones it gives good results. I still don't know why, it requires further debugging.
    
    - **METHOD EXPLANATION** 
        - For each points to be tested, we create two lists from the same starting list of candidates (*struct2*). The first one is sorted by the X coordinates and the second by the Y coordinates.
        - We insert the originating point in the lists at the right places and we take its index.
        - We test a certain number of points around the middle point created by the originating point and put them in two separated lists for X and Y coordinates. This number of points represent the length of the submatrice tested. This "*search_depth*" (aka. the number of points taken around the middle point) is determined by the length of the array tested (length <10 : depth = length/2 ELSE depth=10). 
        - We then keep the intersection between these two arrays and now we have three use cases:
            1. There are no points in the intersection array: we start again with a bigger submatrice.
            2. There is a single point in the intersection : this is the closest point and we keep it.
            3. There are several points : we test the individual distances and keep the smallest one.
 
    - We found out why the results given by bigger matrices are wrong.
    
        - With a simple visualisation tool, we highlighted the points where the distance calculated is wrong, i.e. where the search for the closest point is wrong:
        ![Example](STRUCT1_+_STRUCT2_+_DIFF.png)
        
        - In this screenshot, we can see:
            - In red and blue : the points respectively from the first structure and the one compared.
            - in light blue and light red are the distances correctly calculated.
            - In green : the common points between the two structures.
            - The arrows point to the nearest point calculated by our new method.
            
        - We now give a new screenshot which is centered around the point [413,457]
            ![Example](X_-_Y_SEARCH_EXAMPLE_FOR_A_WRONG_POINT.png)
            - In green : the search list sorted by the X coordinates.
            - In yellow : the search list sorted by the Y coordinates.
            - The intersection is pointed to by the arrow from the originating point.
                - Unfortunately, it is not the nearest point. 
                - This is due to the immediate proximity of a cluster of points which flaws the search.
            - To solve this problem, we need a bigger search depth. (It was set to ten)
        - for the example above, a search depth of a minimum of 26 is required. (tested manually)
        - How to determine a good value for the search depth?
            - We have some ideas : 
                - Since clusters of points seem to indicate for a larger depth, we could find a way to calculate an index for points disparity.
                - This index needs to be calculated quickly and reliably.
    
        - For now, we don't know how to calculate this depth.
            
        - **Successfuly implemented double list search.**
        
        - The structures used in this example are:
        "(((((......(...[..).......)))).).........((....))...((((......((.]...)).....))))...((((.........(((((..(((.....)))..))))....)........(((((..(((...(((.((((((((((....({((((......)))))[[(...).(((....}]]))).......)))))))..))))))...)))..)))).)...)....(((((((........))))))).....(((((........)))))..)))..................................((....))............(..((((....))))..........).....................(.(((((((....((...((((.(((((.....(((....(((((.(((.((((...((((((.(......)))))))..)).)))..)))))))....)))....))))).))))..((((((((.(((((((.((((((......(((..((((..(((....))).))))...)))....)))))).))))))).....(((....))).))))))))...))..))))))).)..((((........))))..........(((.(((((((((..(((...((((......))))...)))...((....))........(((((...(((....)))....))))).(((((((....(((...)))....)).))).))....))))))))).)))..((.(...((((((.....((((((((((.(((....(((((((....))))))).....)))...(((..((.((........)))).)))..).))))))))).(.((.......))...)..))))))......).))............(.(((....))).)......." 
        and
        "((((((...................))))).).........((....))...((((......((.....)).....))))...((((.........((..(((((.....)))))..).)............((((((...((.((..((((((((.(....((((........))))[[(...).(((.....]])))....).)))))..)))..)).....))..))))))....)....((((..............)))).....(((((......)))))..)))..................................((....))...............((((....))))..................................((((....((.......((.........(..(((((.....))))).)..............))......(..(((.(((.((((..((((.........(((((....(((((((((..)))))))..))...))))).......))))...)))).))).).......(((....))).....))..)........))....))))...(((........)))..............(((.(((..(..(.(...((.....))...).)..................(.(((...((((....))))....))).)...(.((((....(....)....))).).).......)..))).))).....((.(...((((((.....(((((((((....((...((((.(......).))))...).).....(.(..((((........)))).).)....)))))))))....................))))))......).))...........(.(.(....).).).........................."
        
        They are purposefully complicated and long to represent a worst case scenario.
        
- The new verbose mode:  
    - Speed of the algorithm with verbose mode enabled:
        - Since the old approach printed out the nearest point found as it was calculated, the calculation needed to wait for the sentences to be printed (this meant that a 0.5s calculation would take ~6s when everything was printed out)
        - With the new approach, the nearest point found are printed at the end of the calculation with the distance value.
        - This allows for a more efficient calculation in verbose mode and it nearly doesn't impact peformances.
    - What is said:
        - When verbose mode is enabled, it indicates when a pool of cores is created and terminated and with how much cores.
        - At the end of the calculation, every nearest point found is printed neatly with an indication of the cache used and the distance calculated.
        - When inverting the structures, it indicates that the algorithm is creating a new pool of cores and changing the base structure.

- Added a random weighted structure file generator.
    - This creates a .fa file with a certain number of structures (default: 100).
    - There are parameters to control the appearance of the structures (density, bias...)
    
- Further optimised the program to create a unique pool of cores for all calculations.


- **ADVICE ON NUMBER OF CORES TO USE**
    
    - When testing a large number of big structures we recommend using a maximum of *6* cores.
    
    - When testing a large number of small structures, we recommend *4* cores.
    
    - When testing a very small number of structures, without taking into account the scale of the structures, always use *2* cores.
    
    - *SPEED COMPARISON*:
        - For a random weighted structure file with 100 structures of length 500:
            -  With 6 cores, we go from 20s of calculation with AptaMat to 6s with AptaFast.
            
        - In general, expect a 2 to 5 times improvement on speed regardless of file size, structure number and length.
        
- **SEARCH_DEPTH CALCULATION**
    
    - We chose to determine empirically the best depth for each matrix size. How did we do it?
        1. We tested the accuracy of AptaFast for matrix sizes of 10, 25, 50, 75, 100, 250, 500 and 1000.
            - For each size, we tested the search_depth limit at which the program started to give improper results.
                - We used a randomly generated set of 1000 dotbracket structures for each size. (density = 0.5; bias = 0.6)
        2. We plotted the result on a graph and chose a linear regression to have a depth for each size value.
            - Here is this plot with all the regression parameters:
            !["Graph"](LINEAR_REGRESSION_GRAPH.png)
            !["Parameters"](LINEAR_REGRESSION_PARAMETERS.png)
            
            - The exact points plotted were:
            
            |X|Y|
            |:-----:|:-----:|
            |10|3|
            |25|4|
            |50|5|
            |75|5|
            |100|6|
            |250|7|
            |500|9|
            |1000|13|
            
        3. We now have a mathematical equation to calculate the search_depth with respect to the length of the structure: *depth = 0.009125 x length + 4.207*
            - For security purposes, we added 2 to the result to counter the rounding error.
        
    - Adding a new `-speed` parameter.
        - When using the program on 1000x1000 non-aligned matrices taken from a real database, we found a minimal search depth to be 26, not 13 as given by the equation.
            - That is why we decided to add a "speed" parameter defaulted to "slow" which determines the greedyness of the algorithm regarding the depth.
            
        - When set to slow (**recomended**) the program doubles the search_depth given by the equation to counter extreme cases.
        
        - Why the name 'speed' and not 'greedyness' ?
            - because lowering the depth of search theoretically improves the calculation speed.
            - In reality it *MIGHT* improve speed by a little when comparing very big structures.
                - When testing a file with 10000 weighted structures of length 100 on 6 cores, we saw no difference in speed whether the speed was set to "slow" or "quick".
                - That is why we chose "slow" as default. 
        
- **WARNING** :
    - It seems important to specify that, for all the performance tests, we used the random structure generator which generates structures impossible to find in real life.

#### FUTURE CHANGES AND IDEAS

- Improving the random structure generation to create realistic dotbracket notation for performance tests.
    
- Taking a look at [GPU optimisation](https://developer.nvidia.com/how-to-cuda-python)

- Taking a look at the placement of gaps and the minimisation of the AptaMat distance.
    - starting development of an alignment tool for the minimisation of the distance.
    
- **IF POSSIBLE**:
    - improving the speed of file reading in the case of very big files.
        - improving ram use in the case of big files.


### WEEK 3 - 16/09/2024 -> 22/09/2024.

- Succesfully divided RAM usage by at least 8 when parsing a file.
    - using uint8 (1byte) instead of default float64 (8bytes) when parsing a file.
    - This greatly improved speed on very big structures (from 214s to 175s for 250 5000-length structures)
    - We see no difference in speed when calculating distance for intermediate to smaller structure for aptafast.
    - Parsing file may be faster but there is a need to find a better way to calculate the dotplot matrix or even completely bypass the matrix and going from dotbracket to coordinates directly.

- Changed the way random_gen.py behaves to generate more realistic structures and to be more reliable and quick.
    - Now, seeing this pattern is impossible : `()`
    - We adjusted the bias and density parameters.
    - the verify function now tests the length of the returned structure.
    
- Supplementary performance test with realistic randomly generated structure files:
    - It seems that with smaller structures, the change in matrix size altered the speed of aptafast in a bad way when comparing it with aptamat which contradicts what we saw last week.
        - You can expect aptafast to be around 1.5 to 2 times slower. We still don't know why it happens and we need to do more tests.
    - For bigger structures, we see a real difference:
        - with 5000 1000-length weighted structure file, we see these results:
        
        |/|SPEED|RAM USAGE|
        |:-----:|:-----:|:-----:|
        |SPEED 'SLOW'|434s|8.2Go|
        |SPEED 'QUICK'|334s|8.2Go|
        |APTAMAT|1077s|45.3Go|
        
    - additionnal tests need to be conducted for smaller structures (under 100-length) and for accuracy.

- Successfuly implemented a faster file parser using a dictionnary and partial search.
    - Expect a 1.33 times improvement of file parsing speed for all files.
    - Accuracy is perfect, expect no coordinates errors.
    
- Completely optimized file parsing with multiprocessing.
    - If "quick" mode is selected, the program WILL use all the cores of the cpu.
    - If "slow" mode is selected, the program will use half the cores of the cpu.
    
    - We saw an immediate improvement in overall speed:
        - for a file with 1000  1000-length structure, we went from 130.21s to 61.37s (quick mode and 6 cores for calculations)
            - The file parsing took 5s instead of 70s.

- **First finished version of AptaFast without gpu optimization.**

**CLEANING**

- Deleted Cache because it slowed things down.

- Cleaned the program for a better efficiency.
    - Now, calculation_core and calculation_core_naive handles the verbose mode and only returns the distance found, not the point.
        - There is now no need to recalculate the distances when the pool of cores has finished.
    - Verbose mode can be very unorganised since all the cores dumps their results on the screen as it is calculated.
    
- Using Naive search because of better performances for structures under 250 of length. (please see performances tests below)

##### NEW PERFORMANCE TESTS:

We tested 67 use cases with files of 10 000 structures of length ranging from 100 to 1000. All the times are in seconds.
We also tested with 2, 4 and 6 cores, every time with the "SLOW", "QUICK" or naive method.

!["With Line Correlation"](new_perf_test_line_corr.png)

- In this screenshot of the tests, we have colored the cells line per line, so colors from any line are completely independant from another.
    - We can see here how and when AptaFast becomes faster than AptaMat in each use cases. And also where using naive search in Aptafast is completely equivalent to the double list search.
    - We can also see when naive search becomes slower: at around 250 in length, this is where we switch methods in the program.
    
    
!["With all correlations"](new_perf_test_all_corr.png)

- In this one, the colors in all cells are corelated. We can see how AptaFast is more homogeneous in its speed than AptaMat.

- How to explain the sudden pics in time ?
    - Since we used the same files for each lines, we can assume that we encountered some kind of "structural extreme" where the program encounters a worst case scenario for one of its mode (hence the pics for 2 cores at 250 length).
    
    
**TEMPORARY CONCLUSION ON PERFORMANCE**

- AptaMat is substantially faster than aptafast for smaller structures.
    - This surely comes from the fact that singlethreaded performance is very important in the case of AptaMat. 
    - Moreover, when calculating with a very small number of points, the naive implementation is still faster because there is a sequence of very small and efficient operations that are very fast when running sequentially.
    - When calculating bigger structures, the double-list method becomes faster because the number of points to be tested in the naive implementation is enormous.
    

- AptaFast's naive implementation is not conclusively faster than the double list search for smaller structures.
    - It is still bound by the management of CPU cores, even if it's the exact same method used in AptaMat.
    - We still want to use it for smaller structure since the double list search becomes completely identical to the naive search (especially in the "SLOW" mode).
        - Why? : Because there is less instructions and less variables manipulated by the program in naive mode.
      
    
#### FUTURE CHANGES AND IDEAS

- **Making a definitive aptafast implementation with:**
    - multiprocessing
    - GPU optimizing using PyCuda
    - Naive + double list search

- **HOW:**
    - Having the search for the nearest point be done by the GPU for each point.
    - Implementing multiprocessing at a higher level:
        - At the structure level and not point level.
            - This means that we would allocate the comparison for a structure on a single thread which will call the GPU for the nearest point search.
    - Using Naive search or double list search with a threshold using the GPU.


### WEEK 4 - 23/09/2024 -> 29/09/2024.

**First Release of AptaFast**

- We have know a definitive implementation of AptaFast with overall algorithmic optimization and multiprocessing.

- Some improvement could still be made regarding multiprocessing depth. Do we use it on a *per point* basis or *per structure* basis ? Testing still needs to be made.

- Accuracy testing needs to be conducted on a more automated and higher level.

**No GPU optimization**
- We decided to abandon GPU optimization since it would be a problem for portability.

**Definitive accuracy and speed testing**

- Here is the updated performance test with an accuracy test:

!["Perf + acc test"](new_perf_acc_tests.png)

- We tested 10 new sets of 100 to 1000 length structures. There were 1000 structures in each of them.

- The time now takes into account the file parsing time as well as the calculation time.

- Accuracy and average deviation are calculated by taking to average of the accuracy (in %) and the average of the average devation for each pass.

- We can see that, the bigger the structure is, the worse the accuracy becomes. BUT, the average deviation stays really small which is satisfactory.

- SO:
    - If you want the best accuracy and speed, we reccomend using the "SLOW" mode with 6 cores.
    - If you want the best speed regardless of perfect accuracy, use the "QUICK" mode with 6 cores.

**Changes in clustering**

- Started understanding, fixing and cleaning clustering.

- Updated file paths to Aptafast and aptamat.

- Added CPU and GPU accelerated functions to view the affinity matrix using matplotlib and vispy.
    - GPU function still work in progress.
    
- Cleaned a bit more.

- Finished GPU accelerated affinity matrix visualizator.
    - We can see the figures created by both histogram methods:
    
!["CPU"](3D_histogram_CPU.png)

- The histogram created with matplotlib, CPU accelerated.
    
!["GPU"](3D_histogram_GPU.png)

- The histogram created with vispy, GPU accelerated with an openGL implementation.
    
- The CPU implementation is very slow especially to turn in 3D space.

- The GPU implementation may be a bit less readable but faster to display and explore.

- Adapted AptaFast to work with the clustering algorithm.
    - Clustering is two times faster.


- Finished cleaning up the clustering algorithm.

- `clustering_AptaMat.py` is to be used in the command line.

#### FUTURE CHANGES AND IDEAS

- Adding a new way to visualize data with the clustering algorithm because the affinity matrix is very big and not very readable.

- Implementing an alignement tool to minimize the AptaMat distance and making AptaMat more independant.
    - First, naive implementation by testing all the possibilities.
    - Then possible algorithmic optimization.
    - In the end, maybe using a genetic algorithm.
