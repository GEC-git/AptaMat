#! /usr/bin/env python3
import argparse
import string
import os, sys
import numpy as np
from scipy.spatial.distance import cityblock
from warnings import warn as Warning
from typing import List


class Parse_File:
    """
    Structure parser.

    Take as input either file containing structures or dotbracket structure as `str`
    """

    def __init__(self, _input: str):
        self.header = None
        self.sequence = None
        lines = open(_input, "r").readlines()

        if not lines:
            raise ValueError("Empty file.")

        self.structures = []
        self._parse(lines)

    def _parse(self, lines: List[str]):
        """
        Parse dotbracket structures from an handled file object
        :param lines:
        :return:
        """

        for i, line in enumerate(lines):

            if line.isspace():
                continue

            if line.startswith('>'):
                self.header = line.strip('\n >')
                continue

            if line.strip().isalpha():
                self.sequence = line.strip()
                continue

            if line.startswith('--'):
                id = line.strip('\n -')
                db = lines[i + 1].strip('\n')
                self.structures.append(Dot_Bracket(db, id))

            else:
                if lines[i - 1].startswith('--'):
                    continue
                db = line.strip('\n')
                self.structures.append(Dot_Bracket(db))

        return self.structures


class Dot_Bracket:
    """ Create a DotBracket object"""

    def __init__(self, dotbracket: str, id: str = None):

        if any('(' == char or '.' == char for char in dotbracket):
            self.db = dotbracket
        else:
            raise ValueError("Bad structure format \nBe careful that the input is in dotbracket format ")

        self.id = id
        self.dotplot = Dotplot(self)

    def __eq__(self, other):
        if isinstance(other, Dot_Bracket):
            return self.db == other.db


class Dotplot:
    """ Create a Dot_plot from a Dotbracket object"""

    def __init__(self, structure: Dot_Bracket):
        self.matrix = self._build_matrix(structure)
        self.coordinates = self._coordinates(self.matrix)

    def __eq__(self, other):
        if isinstance(other, Dotplot):
            return self.matrix == other.matrix and self.coordinates == other.coordinates

    @staticmethod
    def _build_matrix(structure: Dot_Bracket):
        """
        Convert `Dot_Bracket` structure to matrix format structure.

        Work fine with extended dotbracket notation.

        :return: matrix
        """
        opening = '([{<' + string.ascii_uppercase
        closing = ')]}>' + string.ascii_lowercase
        matrix = np.zeros((len(structure.db), len(structure.db)))

        for i, db1 in enumerate(structure.db):
            open_db = 0

            if db1 in opening:

                for j, db2 in enumerate(structure.db):

                    if j < i:
                        continue

                    elif j == '.':
                        continue

                    elif db2 == db1:
                        open_db += 1

                    elif closing.find(db2) == opening.find(db1) and open_db > 0:
                        open_db -= 1

                    if closing.find(db2) == opening.find(db1) and open_db == 0:
                        matrix[i, j] = 1
                        break

        matrix = np.asarray(matrix)

        return matrix

    @staticmethod
    def _coordinates(matrix: np.ndarray):
        """
        Translate binary matrix (dotplot) into coordinate array.
        1 -> dot in the matrix

        :return: coordinates
        """
        coordinates = []
        for r, row in enumerate(matrix):
            # print(r)
            for c, col in enumerate(row):
                if col == 1:
                    # print(c)
                    dot = [r + 1, c + 1]
                    # print(dot)
                    coordinates.append(dot)
        coordinates = np.asarray(coordinates)

        return coordinates


class NucleicAcid:
    Nuc = ('G', 'A', 'T', 'U', 'C')
    Comp = ('C', 'T', 'A', 'A', 'G')

    def __init__(self):
        pass


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-structures',
                        nargs='+',
                        type=str,
                        help='2D structures in dotbracket notation')
    parser.add_argument('-files',
                        action="store",
                        nargs='+',
                        help='Input file containing structures to be compare. The first is considered as the '
                             'Template')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(args.structures) < 2:
        raise ValueError('Missing one argument in -structures')

    return args


def file_process(files: str):
    """
    Applied with file(s) as argument. It iterate over the input files to extract the dotbracket notation and calculate
    the AptaMAT between them.

    Note that it is recommended to either put all the structures in one file, or one structure per file::

        >Example_1
        SEQUENCE (optionnal)
        --Structure--
        .(((............))).


    :param files: input's filepath
    """

    for n, file in enumerate(files):

        parsed_file = Parse_File(file)
        print(parsed_file.header)

        for i, compared_struct in enumerate(parsed_file.structures):

            if i == 0 and n == 0:
                template_struct = compared_struct

            else:
                if template_struct.id is not None or compared_struct.id is not None:
                    print(template_struct.id, '-', compared_struct.id)

                print('', template_struct.db)
                print('', compared_struct.db)
                distance_Manhattan(template_struct, compared_struct)


def distance_Manhattan(struct_1: Dot_Bracket, struct_2: Dot_Bracket):
    """
    Calculate distance between struct_1, struct_2 using Manhattan distance
    with ::
        - struct_1 : Template structure
        - struct_2  : Compared structure

    The applied formula is the following :

    :math:`Distance (struct 1, struct 2) + Distance (struct 2, struct 1)/(n distance (struct 1, struct 2) + n
    distance (struct 2, struct 1) / ( Length(struct 1) - 4 )`


    Where n distance (struct 1, struct 2) is the number of minimum distance handled between struct 1 and struct 2, and
    n distance (struct 2, struct 1) the same with struct 2 and struct 1

    :param struct_1: Dot_Bracket object referring to the template structure
    :param struct_2: Dot_Bracket object referring to the compared structure
    :return: dist: Distance calculated

    """

    if len(struct_1.db) != len(struct_2.db):
        raise Exception("Compared structures must have the same length")

    if struct_2 == struct_1:
        print('> AptaMAT : 0\n')
        return 0

    elif struct_1.dotplot.coordinates.size > 0 and struct_2.dotplot.coordinates.size > 0:

        # Compare Predicted --> Original structure
        dist_PO, nb_compare_PO = compute_distance(struct_1, struct_2)

        # Compare Original --> Predicted structure
        dist_OP, nb_compare_OP = compute_distance(struct_2, struct_1)

        dist = (dist_OP + dist_PO) / (nb_compare_PO + nb_compare_OP)
        print(f'> AptaMAT : {dist}\n')
        return dist

    elif struct_1.dotplot.coordinates.size == 0:
        Warning("Template structure is not folded", stacklevel=2)

        return None

    else:
        Warning("Compared structure is not folded", stacklevel=2)

        return None


def compute_distance(struct_1, struct_2):
    """
    Proceed to distance analysis between input struct_1 and struct_2

    :param struct_1: 2D np.array matrix presenting coordinates
    :param struct_2: 2D np.array matrix presenting coordinates
    :return: distance, compared_point
    """
    nearest_dist = []

    for point_1 in struct_1.dotplot.coordinates:
        point_dist = []
        for point_2 in struct_2.dotplot.coordinates:
            # print(point_P, point_O)
            Man_Dist = cityblock(point_1, point_2)
            point_dist.append(Man_Dist)

        if point_dist:
            nearest_dist.append(min(point_dist))

    distance = sum(nearest_dist)
    compared_point = len(nearest_dist)

    return distance, compared_point


if __name__ == '__main__':
    args = parse_arguments()

    if 'structures' in vars(args).keys():
        for i, structure in enumerate(args.structures):
            if i == 0:
                _structure_t = Dot_Bracket(structure)
            else:
                _structure_c = Dot_Bracket(structure)
                print('', _structure_t.db)
                print('', _structure_c.db)
                distance_Manhattan(_structure_t, _structure_c)

    elif 'files' in vars(args).keys():
        file_process(**vars(args))

    else:
        sys.exit(1)
