###############################################################
#           APTAFAST FOR RECURRING USE OF THE CALCULATIONS
###############################################################


#! /usr/bin/env python3
import argparse
import os
import pathlib
#import string
import sys
import warnings
from copy import deepcopy
import time as tm
import multiprocessing as mp
import numpy as np
from scipy.spatial.distance import cityblock, euclidean

np.set_printoptions(threshold=sys.maxsize)


#############################################################
#           CLASS
#############################################################

class Parse:
    """
    Structure parser from fasta or alignment formatted files
    """
    suffix = ('.aln', '.a', '.fa', '.fasta', '.db', '.txt')

    def __init__(self, file):

        if pathlib.Path(file).is_file():
            self.path = os.path.abspath(file)
        else:
            raise FileNotFoundError

    @classmethod
    def file(cls, file, pool:object ,verbose=False):
        """
        Check file extension before parsing
        """
        if pathlib.Path(file).suffix in cls.suffix:
            if verbose:
                print(f'Reading file: {file} ...\n')
            out = cls.infile(file, pool, verbose)
        else:
            raise AttributeError('{} extension not compatible. Please try with {}\n'.format(pathlib.Path(file).suffix,
                                                                                            ', '.join(cls.suffix)))
        return out

    @classmethod
    def infile(cls, file, pool:object, verbose=False):
        """
        Structure parser from fasta formatted files

        Examples
        --------

        >ID\n
        TCGATTGGATTGTGCCGGAAGTGCTGGCTCGA\n
        --Template--
        ((((.........(((((.....)))))))))\n
        [ weight ]\n
        --Compared1--\n
        .........(((.(((((.....))))).)))\n
        [ weight ]--\n
        --Compared2\n
        ..........((.((((.......)))).)).\n
        [ weight ]\n

        Parameters
        ----------
            file: pathlib.Path
                File to be parsed.
            verbose:
                True or False.

        Returns
        -------
        :return: structures list
        """

        if pathlib.Path(file).suffix not in cls.suffix:
            raise AttributeError(
                'Tried to parse a fasta formatted file but got {} instead.\n'.format(pathlib.Path(file).suffix))

        lines = open(file).readlines()
        structures = []
        non_parsed_struct=[]
        non_parsed_weights=[]
        for i, line in enumerate(lines):
            if verbose:
                print(line)

            if line.isspace():
                continue

            if line.startswith(('>', '--')):
                id = line.strip('\n >-')
                if not Dotbracket.is_dotbracket(lines[i + 1]):
                    sequence = lines[i + 1].strip('\n')
                else:
                    sequence = None
                continue

            elif Dotbracket.is_dotbracket(line):
                dotbracket = line.strip()
                non_parsed_struct.append((dotbracket, sequence, id, file))
                #struct = SecondaryStructure(dotbracket, sequence, id, file=file)
                try:
                    lines[i + 1]
                except IndexError:
                    pass
                else:
                    if lines[i + 1].startswith(('>', '--')):
                        pass
                    elif not Dotbracket.is_dotbracket(line):
                        non_parsed_weights.append(float(line))
                        #struct.weight = float(lines[i + 1].split(sep=' ')[1])
                    else:
                        non_parsed_weights.append(0)
                        #struct.weight = 0

            else:
                continue
           
        structures=pool.starmap(SecondaryStructure, [(non_parsed_struct[i][0], non_parsed_struct[i][1], non_parsed_struct[i][2], non_parsed_struct[i][3]) for i in range(len(non_parsed_struct))])
            
        def get_id(struct):
            num=""
            for char in struct.id:
                if 48<=ord(char)<=57:
                    num+=char
            return int(num)
        
        structures=sorted(structures, key=lambda struct : get_id(struct))
        
        for i,elt in enumerate(structures):
            elt.weight = non_parsed_weights[i]
                
        #structures.append(struct)

        if not structures:
            return AttributeError("File content not compatible (malformed, corrupted or empty\n"
                                  "Check documentation to know which type of file is supported\n")

        return structures

class Dotbracket:
    """Create a DotBracket object"""
    #gap_penalty_matrix = [0,0]
    gap_penalty_matrix = [1,1]
    def __init__(self, dotbracket: str, gap_penalty=None):
        self.dotbracket = None
        self.gap = None

        if gap_penalty is None:
            gap_penalty = self.gap_penalty_matrix

        if self.is_dotbracket(dotbracket):
            self.set_dotbracket(dotbracket, gap_penalty)
        else:
            raise ValueError("Bad structure format \nBe careful that the input is in dotbracket format\n")

    def __str__(self):
        return self.dotbracket

    def __eq__(self, dotbracket):
        if isinstance(dotbracket, Dotbracket):
            return self.dotbracket == dotbracket.dotbracket

    def __len__(self):
        return len(self.dotbracket)

    def set_dotbracket(self, dotbracket, gap_penalty=None):
        if gap_penalty is None:
            gap_penalty = self.gap_penalty_matrix
        self.dotbracket = dotbracket
        self.gap = self.gap_penalty(gap_penalty)

    def gap_penalty(self, gap_penalty):
        if not gap_penalty[0]:
            return 0
        else:
            return self.dotbracket.count("-")*gap_penalty[0]
        # for i, c in enumerate(self.dotbracket):
        #     if c == '-' and self.dotbracket[i - 1] == '-':
        #         penalty += gap_penalty[1]
        #     elif c == '-':
        #         penalty += gap_penalty[0]
        #     else:
        #         pass


    @staticmethod
    def is_dotbracket(dotbracket):
        if any(char.isdigit() for char in dotbracket):
            return False
        if any('(' == char or '.' == char for char in dotbracket):
            return True


class Dotplot(Dotbracket):
    """ Create a Dotplot from a dotbracket object"""

    def __init__(self, dotbracket):
        Dotbracket.__init__(self, dotbracket)
        self.matrix = self.build_matrix()
        self.coordinates = self.get_coordinates()

    def __str__(self):
        return self.matrix

    def __eq__(self, other):
        if isinstance(other, Dotplot):
            return self.matrix == other.matrix and self.coordinates == other.coordinates

    def __len__(self):
        return self.matrix.size, self.matrix.shape

    def build_matrix(self):
        """
        Convert `Dotbracket` structure to matrix format structure.

        Works fine with extended dotbracket notation.

        Returns
        --------
            matrix: np.ndarray
                binary matrix.
        """
        # opening = '([{<' + string.ascii_uppercase
        # closing = ')]}>' + string.ascii_lowercase
        matrix = np.zeros((len(self.dotbracket), len(self.dotbracket)),dtype=np.uint8)
        par="()"
        cro="[]"
        acc="{}"
        #void=".-"
        dict_par={}
        dict_cro={}
        dict_acc={}
        for i,elt in enumerate(self.dotbracket):
            if elt in par:
                dict_par[i]=elt
            elif elt in cro:
                dict_cro[i]=elt
            elif elt in acc:
                dict_acc[i]=elt
        
        while dict_par != {}:
            index_list=list(dict_par.keys())
            i=0
            while dict_par[index_list[i]]==dict_par[index_list[i+1]]:
                i+=1
            matrix[index_list[i],index_list[i+1]]=np.uint8(1)
            del(dict_par[index_list[i]])
            del(dict_par[index_list[i+1]])
        try:    
            while dict_cro != {}:
                index_list=list(dict_cro.keys())
                i=0
                while dict_cro[index_list[i]]==dict_cro[index_list[i+1]]:
                    i+=1
                matrix[index_list[i],index_list[i+1]]=np.uint8(1)
                del(dict_cro[index_list[i]])
                del(dict_cro[index_list[i+1]])
                
            while dict_acc != {}:
                index_list=list(dict_acc.keys())
                i=0
                while dict_acc[index_list[i]]==dict_acc[index_list[i+1]]:
                    i+=1
                matrix[index_list[i],index_list[i+1]]=np.uint8(1)
                del(dict_acc[index_list[i]])
                del(dict_acc[index_list[i+1]])
        except IndexError:
            #In the case of a non closing pseudoknot.
            dict_cro={}
            dict_acc={}
        return matrix

    def get_coordinates(self):
        """
        Translate binary matrix (dotplot) into coordinate array.
        1 -> dot in the matrix

        Returns
        --------
            coordinates: np.ndarray
                array of coordinates (x,y)
        """
        coordinates = []
        for r, row in enumerate(self.matrix):
            # print(r)
            for c, col in enumerate(row):
                if col:
                    dot = [r + 1, c + 1]
                    coordinates.append(dot)
        coordinates = np.asarray(coordinates)
        return coordinates


class SecondaryStructure(Dotplot):
    """
    Create a Secondary Structure object from a dotbracket (str)

    Secondary Structure object shares its instances with Dotbracket and Dotplot classes
    """

    def __init__(self, dotbracket: str, sequence: str = None, id: str = None, file=None, fam=None):
        self.id = id
        self.sequence = sequence
        self.family=fam
        Dotplot.__init__(self, dotbracket)
        self.weight = 0
        if file is not None:
            if not os.path.isfile(os.path.abspath(file)):
                raise AttributeError("File {} does not exist.\n".format(file))
            else:
                self.file = file
                
    def __eq__(self, ss):
        if isinstance(ss, SecondaryStructure):
            return self.dotbracket == ss.dotbracket

    def __len__(self):
        return len(self.dotbracket)


#############################################################
#           Internal Functions
#############################################################

def _warning(msg, *args, **kwargs):
    return str(msg) + '\n'


warnings.formatwarning = _warning


def _check_instance(obj):
    if isinstance(obj, SecondaryStructure):
        pass

    elif isinstance(obj, Dotplot):
        pass

    elif isinstance(obj, Dotbracket):
        obj = Dotplot(obj.dotbracket)

    elif isinstance(obj, str):
        obj = SecondaryStructure(obj)

    else:
        raise AttributeError('Object of type {} not recognized'.format(type(obj)))

    return obj


def _result_print(template_struct, compared_struct, weight=None):
    print(template_struct.id, end=' ')
    print('-', compared_struct.id)
    if weight is not None:
        print(f'Weigth= {str(weight)}')
    print('', template_struct.dotbracket)
    print('', compared_struct.dotbracket)
    print('> AptaMat:\n', end='  ')


#############################################################
#           Functions
#############################################################

def compute_distance_clustering(struct_1: object, struct_2: object, method, speed:str, verbose=False):
    r"""
    Calculate distance between struct_1, struct_2 using Manhattan distance without using multiprocessing as it is used at a higher level in clustering_AptaMat.
    with ::
        - struct_1 : Template structure
        - struct_2 : Compared structure

    The applied formula is the following

    .. math::
        D_AM(struct_1, struct_2) = \frac{\sum_{P\in P_A}d(P,P_B) + \sum_{P\in P_B}d(P,P_A) + N_G}{Card(P_A) + Card(P_B)}

    Where, for any given point :math:`P = (x, y) ∈ R^2` and any finite subset
    :math:`C ⊂ R^2`, we denote by :math:`Card(C)` the cardinal of C, and by :math:`d(P,C)` the
    Manhattan distance from P to its nearest neighbor in C.
    Where n distance (struct 1, struct 2) is the number of minimum distance handled between struct 1 and struct 2, and
    n distance (struct 2, struct 1) the same with struct 2 and struct 1

    Parameters
    ----------
    struct_1 :  SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    struct_2: SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    method: str
        Method for distance calculation.
    verbose :
        True or False
    cache :
        CompressedCache passed through

    Returns
    -------
        dist : float or None
            The AptaMat distance between struct_1 and struct_2. None if structures are not folded
    """

    # Check input type before running calculation
    s1 = deepcopy(struct_1)
    s1 = _check_instance(s1)
    s2 = deepcopy(struct_2)
    s2 = _check_instance(s2)

    if verbose:
        print(f'Comparing :\n'
              f'{struct_1.id}\n'
              f'{struct_1.dotbracket}\n'
              f'{struct_2.id}\n'
              f'{struct_2.dotbracket}\n')

    # Check length of compared structures
    # Alignment is strongly recommended
    
        if len(s1) != len(s2):
            warnings.warn(
                "Input structures with different sizes.\n "
                "For accurate results, please perform sequence or structure alignment before. \n")

    if s2 == s1:
        return 0

    elif s1.coordinates.size > 0 and s2.coordinates.size > 0:

        # Template structure --> Compared structure
        nearest_dist=[]
        
        struct2=[list(elt) for elt in struct_2.coordinates]
        struct1=[list(elt) for elt in struct_1.coordinates]
        if len(struct_2) >=200:
            search_depth=30
            for elt in struct1:
                nearest_dist.append(calculation_core(elt, struct2, method, search_depth, verbose))
        else:
            for elt in struct1:
                nearest_dist.append(calculation_core_naive(elt, struct2, method, verbose))
        
        d_TC = sum(nearest_dist)
        
        nearest_dist=[]
        
        # Compared structure --> Template structure
        if len(struct_1) >=200:
            search_depth=30
            for elt in struct2:
                nearest_dist.append(calculation_core(elt, struct1, method, search_depth, verbose))
        else:
            for elt in struct2:
                nearest_dist.append(calculation_core_naive(elt, struct1, method, verbose))
        
        d_CT = sum(nearest_dist)

        dist = (((d_CT + d_TC) + (s1.gap + s2.gap)) / (len(s1.coordinates) + len(s2.coordinates)))
        if verbose:
            print(f'({d_CT} + {d_TC}) + ({s1.gap} + {s2.gap})) / ({len(s1.coordinates)} + {len(s2.coordinates)})')
        return dist

    elif s1.coordinates.size == 0:
        warnings.warn("Template structure is not folded."+str(struct_1.id)+str(struct_1.dotbracket))
        return None

    else:
        warnings.warn("Compared structure is not folded."+str(struct_2.id)+str(struct_2.dotbracket))
        return None


def compute_distance(struct_1: object, struct_2: object, method, nb_pool: int, pool: object, speed:str, verbose=False):
    r"""
    Calculate distance between struct_1, struct_2 using Manhattan distance
    with ::
        - struct_1 : Template structure
        - struct_2 : Compared structure

    The applied formula is the following

    .. math::
        D_AM(struct_1, struct_2) = \frac{\sum_{P\in P_A}d(P,P_B) + \sum_{P\in P_B}d(P,P_A) + N_G}{Card(P_A) + Card(P_B)}

    Where, for any given point :math:`P = (x, y) ∈ R^2` and any finite subset
    :math:`C ⊂ R^2`, we denote by :math:`Card(C)` the cardinal of C, and by :math:`d(P,C)` the
    Manhattan distance from P to its nearest neighbor in C.
    Where n distance (struct 1, struct 2) is the number of minimum distance handled between struct 1 and struct 2, and
    n distance (struct 2, struct 1) the same with struct 2 and struct 1

    Parameters
    ----------
    struct_1 :  SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    struct_2: SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    method: str
        Method for distance calculation.
    verbose :
        True or False
    cache :
        CompressedCache passed through

    Returns
    -------
        dist : float or None
            The AptaMat distance between struct_1 and struct_2. None if structures are not folded
    """

    # Check input type before running calculation
    s1 = deepcopy(struct_1)
    s1 = _check_instance(s1)
    s2 = deepcopy(struct_2)
    s2 = _check_instance(s2)

    if verbose:
        print(f'Comparing :\n'
              f'{struct_1.id}\n'
              f'{struct_1.dotbracket}\n'
              f'{struct_2.id}\n'
              f'{struct_2.dotbracket}\n')

    # Check length of compared structures
    # Alignment is strongly recommended
    if len(s1) != len(s2):
        warnings.warn(
            "Input structures with different sizes.\n "
            "For accurate results, please perform sequence or structure alignment before. \n")

    if s2 == s1:
        return 0

    elif s1.coordinates.size > 0 and s2.coordinates.size > 0:
        if verbose:
            print("Template structure --> Compared structure")
        # Template structure --> Compared structure
        
        if nb_pool>mp.cpu_count() or nb_pool <=0:
            return ValueError("Incorrect number of cores")
        d_TC = pairwise_distance_optimised(s1, s2, method, pool, speed, verbose)
        if verbose:
            print("Compared structure --> Template structure\n")
        # Compared structure --> Template structure
        d_CT = pairwise_distance_optimised(s2, s1, method, pool, speed, verbose)

        dist = (((d_CT + d_TC) + (s1.gap + s2.gap)) / (len(s1.coordinates) + len(s2.coordinates)))
        if verbose:
            print(f'({d_CT} + {d_TC}) + ({s1.gap} + {s2.gap})) / ({len(s1.coordinates)} + {len(s2.coordinates)})')
        return dist

    elif s1.coordinates.size == 0:
        warnings.warn("Template structure is not folded.\n")
        return None

    else:
        warnings.warn("Compared structure is not folded.\n")
        return None

def calculation_core(point1, struct2, method, search_depth, verbose):
    """
        Independant function for determining the closest point
    
        This function uses double list search - please read changelog for further explanations.
        
        Returns a tuple with:
            - The point originating the search
            - The closest point found to the first one
    """

    def Y(x): return x[1]
    def X(y): return y[0]
    
    X_axis=struct2[:]
    Y_axis=struct2[:]
    
    struct2_str=[str(elt) for elt in struct2]
    
    if str(point1) in struct2_str:
        if verbose:
            print("Found nearest point",point1," - ",point1,"| Value:",0)
        return 0

    else:
        search_num=0
        X_axis.append(point1)
        Y_axis.append(point1)
        X_axis = sorted(X_axis, key=lambda pt : X(pt))
        Y_axis = sorted(Y_axis, key=lambda pt : Y(pt))
        X_axis_str=[str(elt) for elt in X_axis]
        Y_axis_str=[str(elt) for elt in Y_axis]
        
        
        X_index=int(X_axis_str.index(str(point1)))
        Y_index=int(Y_axis_str.index(str(point1)))
        
        finished = False
        
        while not finished:
            
            intersection=[]
            
            
            searching_list_X_ascend=X_axis[int(X_index+1):int(X_index+1+search_depth*search_num+search_depth)]
            
            X_descent_index_low=int(X_index-1-search_depth*search_num-search_depth)
            if X_descent_index_low<0 : X_descent_index_low = 0
            
            X_descent_index_high=X_index
            if X_descent_index_high<0 : X_descent_index_high = 0
            
            searching_list_X_descent=X_axis[X_descent_index_low:X_descent_index_high]
            
            searching_list_X = searching_list_X_descent+searching_list_X_ascend
            
            
            searching_list_Y_ascend=Y_axis[int(Y_index+1):int(Y_index+1+search_depth*search_num+search_depth)]
            
            Y_descent_index_low=int(Y_index-1-search_depth*search_num-search_depth)
            if Y_descent_index_low<0 : Y_descent_index_low = 0
            
            Y_descent_index_high=Y_index
            if Y_descent_index_high<0 : Y_descent_index_high = 0
            
            searching_list_Y_descent=Y_axis[Y_descent_index_low:Y_descent_index_high]
           
            searching_list_Y = searching_list_Y_descent+searching_list_Y_ascend
            
            searching_list_Y_str = [str(elt) for elt in searching_list_Y]
            
            
            intersection = [elt for elt in searching_list_X if str(elt) in searching_list_Y_str]
            
            if intersection == []:
                finished = False
                search_num+=1
            else:
                finished=True
                
                
        if len(intersection)==1:
            if method=="cityblock":
                if verbose:
                    print("Found nearest point",point1," - ",intersection[0],"| Value:",cityblock(point1,intersection[0]))
                return cityblock(point1,intersection[0])
            if method=="euclidean":
                if verbose:
                    print("Found nearest point",point1," - ",intersection[0],"| Value:",euclidean(point1,intersection[0]))
                return euclidean(point1,intersection[0])
        else:
            dist_dict={}
            for elt in intersection:
                if method=="cityblock":
                    dist_dict[cityblock(point1,elt)]=elt
                if method=="euclidean":
                    dist_dict[euclidean(point1,elt)]=elt
            
            keep=min(dist_dict.keys())
            if verbose:
                print("Found nearest point",point1," - ",dist_dict[keep],"| Value:",keep)
            return keep

def calculation_core_naive(point1, struct2, method, verbose):
    struct2_str=[str(elt) for elt in struct2]
    if str(point1) in struct2_str:
        if verbose:
            print("Found nearest point",point1," - ",point1,"| Value:",0)
        return 0
    else:
        point_dist = {}
        for point2 in struct2:
            if method == "cityblock":
                Pair_Dist = cityblock(point1, point2)
            if method == "euclidean":
                Pair_Dist = euclidean(point1, point2)
            point_dist[Pair_Dist]=point2
        keep=min(point_dist.keys())
        if verbose:
            print("Found nearest point",point1," - ",point_dist[keep],"| Value:",keep)
        return keep

def pairwise_distance_optimised(struct_1: object, struct_2: object, method, pool: object, speed:str, verbose=False):
    """
    Proceed to the point distance parsing between input struct_1 and struct_2 using
    Manhattan distance.

    The function returns the sum of the nearest distances found for each points.
    
    This function uses double list search, cache and multiprocessing to optimize the algorithm.

    Parameters
    ----------
    struct_1 : SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    struct_2: SecondaryStructure_or_Dotplot_or_Dotbracket
        Input SecondaryStructure, Dotplot or Dotbracket object.
    method: str
        Method for distance calculation.
    verbose :
        True or False
    cpu_cores: integer
        Number of cores chosen. Used for creating the pool of processes

    Returns
    -------
        distance : float
            manhattan distance
    """
    nearest_points=[]
    
    struct2=[list(elt) for elt in struct_2.coordinates]
    struct1=[list(elt) for elt in struct_1.coordinates]
    if len(struct_2) >=250:
        if speed=="quick":
            search_depth=int(0.009125*len(struct_2)+4.207)+2
        elif speed=="slow":
            search_depth=(int(0.009125*len(struct_2)+4.207)+2)*2
        else:
            return ValueError("The speed value is incorrect")
            
        nearest_points.append(pool.starmap(calculation_core, [(struct1[i],struct2,method,search_depth,verbose) for i in range(len(struct1))]))
    else:
        nearest_points.append(pool.starmap(calculation_core_naive, [(struct1[i],struct2,method,verbose) for i in range(len(struct1))]))
    
    nearest_dist=nearest_points[0]
    
    if verbose:
        print("Finished this pass.")
        
    distance = sum(nearest_dist)
    return distance

def adjust_weight(structures, weight):
    """
    Function to adjust weight values in a set of structures.

    ---
    Parameters
    ----------
    structures : list
        List of SecondaryStructure, Dotplot or Dotbracket object.
    weight : list
        List of weights.

    Returns
    -------
        weight : list
            List of adjusted weight
    """
    # warnings.warn('Included weights does not result in ensemble distribution of 100% \n Adjusting weight ...\n')
    if len(structures) > len(weight):
        unfilled = len(structures) - len(weight)

        if sum(weight) < 1:

            proportion = (1 - sum(weight)) / unfilled
            for i in range(len(weight), len(structures)):
                weight.append(proportion)
        else:
            return ValueError(
                "Weights have been included in the set but cannot be completed for {} last structures\n"
                "Check your weight ratio.\n".format(unfilled))

    if len(structures) < len(weight):
        weight = weight[0, len(structures)]

    if sum(weight) > 1:
        weight = [i * 100 / sum(weight) for i in weight]

    return weight

def main():
    parser = argparse.ArgumentParser(description="AptaMat is a simple script which aims to measure differences between "
                                                 "DNA or RNA secondary structures. The method is based on the "
                                                 "comparison of binary matrices representing the secondary "
                                                 "structures.AptaMat algorithm is compatible with the extended\n"
                                                 "dotbracket notation.\n"
                                                 "To compute the calculation, user can input the different dotbracket"
                                                 "structures to be compared, either in raw arguments or input file "
                                                 "argument.\n"
                                                 "The default computation mode is pairwise. Using -ensemble argument, "
                                                 "user can calculate AptaMat distance for an ensemble of structure.")

    parser.add_argument('-v',
                        '--verbose',
                        help="Increase output verbosity.",
                        action="store_true")
    
    parser.add_argument('-speed', 
                        default='slow',
                        help="Using greedy or non greedy depth calculation",
                        nargs='?',
                        choices=['slow','quick'])

    parser.add_argument('-structures',
                        nargs='+',
                        type=str,
                        help='2D structures in dotbracket notation.')
    parser.add_argument('-weights',
                        nargs='+',
                        type=float,
                        help='Additional weight for structures.')
    parser.add_argument('-files',
                        action="store",
                        nargs='+',
                        help='Input files containing structures to be compare. The first is considered as the '
                             'Template.')
    parser.add_argument('-ensemble',
                        action='store_true',
                        help='Calculate AptaMat value for ensemble.')
    parser.add_argument('-method',
                        default='cityblock',
                        nargs='?',
                        choices=['cityblock', 'euclidean'])
    
    
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.structures is not None and args.files is not None:
        raise KeyError('-structures and -files cannot be used at the same time.')


    if args.speed == "quick":
        pooling=mp.Pool(mp.cpu_count())
        print("MODE: quick | using all cores for file parsing")
    else:
        pooling=mp.Pool(int(mp.cpu_count()/2))
        print("MODE: slow | using half the cores for file parsing")
    
    struct_list = []
    weights = []
    file_time_start=tm.time()
    
    ##################################
    #  Input structures preparation  #
    ##################################
    struct_sizes=[]
    if args.structures is not None:
        if len(args.structures) < 2:
            raise ValueError('Missing one argument in -structures')
        if args.weights:
            weights = adjust_weight(args.structures, args.weights)
            
        for i, structure in enumerate(args.structures):
            struct = SecondaryStructure(dotbracket=structure, id='structure' + str(i))
            struct_sizes.append(len(struct))
            if weights:
                struct.weight = weights[i]
            struct_list.append(struct)

    ############################
    #  Input file preparation  #
    ############################
    if args.files is not None:
        if args.weights is not None:
            warnings.warn('Arguments "file" and "weights" are not compatible.\n'
                          'To perform weighted calculation with input file please'
                          'include the values in {}.\n'.format(args.file))

        for file in args.files:
            structures = Parse.file(file,pooling)
            if not sum(struct.weight for struct in structures) == 1:
                weights = adjust_weight(structures,
                                        [struct.weight for struct in structures])
                for struct, w in zip(structures, weights):
                    struct.weight = w
            struct_list += structures

    # Stop the program whether structures input are not valid or not found.
    if not struct_list:
        raise ValueError('No valid structure parsed.\n')
    file_time_finish=tm.time()
    pooling.terminate()
    ##########################
    #  Distance calculation  #
    ##########################

    print("This is a multiprocessed program, you have",mp.cpu_count(),"cores in your CPU.")
    nb=int(input("How much do you want to use? "))
    print("Creating pool on",nb,"cores.\n")
    print("Working...\n")
    start=tm.time()
    pooling=mp.Pool(nb)
    res=[]
    for i, compared_struct in enumerate(struct_list):
        template_struct = struct_list[0]

        if i == 0:
            template_struct.distance = 0

        else:
            compared_struct.distance = compute_distance(struct_1=template_struct,
                                                        struct_2=compared_struct,
                                                        method=args.method,
                                                        nb_pool=nb,
                                                        pool=pooling,
                                                        speed=args.speed,
                                                        verbose=args.verbose)
            if not args.ensemble:
                _result_print(template_struct, compared_struct)
                print(compared_struct.distance, end='\n\n')
                res.append(compared_struct.distance)
    finish=tm.time()
    pooling.terminate()
    if args.verbose:
        print("File parsing time:",round(file_time_finish-file_time_start,2),"s")
        print("Calculation time:",round(finish-start,2),"s")
        tot=round(finish-start,2)+round(file_time_finish-file_time_start,2)
        print("Total: ",tot,"s")
    ##########################
    #  Ensemble calculation  #
    ##########################
    

    if args.ensemble and len(struct_list) > 2:
        set_distance = 0
        for compared_struct in struct_list[1:]:
            if all(s.weight == 0 for s in struct_list):
                if args.verbose:
                    weight = 1 / (len(struct_list) - 1)
                    _result_print(template_struct, compared_struct, weight)
                    print(compared_struct.distance * weight, end='\n\n')
                set_distance += compared_struct.distance * (1 / (len(struct_list) - 1))
            else:
                if args.verbose:
                    _result_print(template_struct, compared_struct, compared_struct.weight)
                    print(compared_struct.distance * compared_struct.weight, end='\n\n')
                set_distance += compared_struct.distance * compared_struct.weight

        print('> AptaMat of structure set \n ', set_distance)

    elif args.ensemble and len(struct_list) <= 2:
        raise ValueError('Not enough input structure to work as a set.\n')
    else:
        pass

if __name__ == '__main__':
    main()
