import sys
import math
import numpy as np


# Node format: (code number, generation, color, east, north, west, south)
# NOTE: sys.maxsize used to denote presence of wall as neighbor
# NOTE: 0 = bottom left, 1 = bottom right, 2 = top left, 3 = top right

def qlao(code, tx, ty, direction):
    """
    This is the circled plus operator from the paper. QLAO stands for quad location addition operator
    :rtype: int
    :param code:
    :param tx: the tx value fom the paper
    :param ty: the ty value from the paper
    :param direction:
    :return:
    """
    return (((code | ty) + (direction & tx)) & tx) | (((code | tx) + (direction & ty)) & ty)


def update_neighbors(tree, cell, problem_size):
    """
    updates of all neighbors of a given cell code regardless of brotherhood. the first for loop of figure 7
    :param tree: linear tree containing neighbors
    :param cell: tuple representing gray cell about to be divided
    :param problem_size: problem size used to generate direction codes
    :return: this is a non-returning functions
    """
    # generate direction increments for finding neighbors
    # NOTE: directions in list are in same order as in tuple, this is used in the loop
    directions = [int('01', 2), int('10', 2), int('01' * problem_size, 2),
                  int('10' * problem_size, 2)]  # west, south, north, east
    neighbor_codes = []
    directions = map(lambda x: x << (2 * (problem_size - cell[1])), directions)
    for i in range(4):
        if cell[3 + i] == sys.maxsize:  # if neighbor is known wall, make no changes to any LD
            # print 'Neighbor {0:d} is wall, doing nothing'.format(i)
            neighbor_codes.append(sys.maxsize)
        else:
            neighbor_codes.append(qlao(cell[0], int('01' * problem_size, 2), int('10' * problem_size, 2),
                                       directions[i]))  # this is theorem 1 from paper
            # print 'Neighbor {0:d} is {1:s}'.format(i, binary_to_quaternary_string(neighbor_codes[i], problem_size))
            neighbor = next((x for x in tree if x[0] == neighbor_codes[i]), None)
            if neighbor is not None and neighbor[1] == cell[1]:
                # print 'Incrementing index %d of %s' % ((3+((2+i)%4)), binary_to_quaternary_string(neighbor[0],
                # problem_size))
                neighbor[3 + ((2 + i) % 4)] += 1

                # print 'Neighbors of {0:s} are : E = {1:s} N = {2:s} W = {3:s} S = {4:s}'.format(
                # binary_to_quaternary_string(cell[0], problem_size), *map(lambda x: binary_to_quaternary_string(x,
                # problem_size), neighbor_codes))


def grey_nodes_in_tree(tree):
    """
    check if there are gray nodes in the tree using built-in next function
    :param tree:
    :return: whether or not there's a gray node in the tree
    """
    return next((node for node in tree if node[2] == 'G'), None) is not None


def binary_to_quaternary_string(code, problem_size):
    """
    convert the int location code to # printable quaternary value
    :param code: binary form of code to be # printed
    :param problem_size: amount of quaternary values
    :return: string form of quaternary
    """
    if code == sys.maxsize:
        return '#'
    desired = ''
    for i in range(1, 2 * problem_size, 2):
        desired += str(int('{0:06b}'.format(code)[i - 1:i + 1], 2))

    return desired


def fetch_color(code, gen, grid, problem_size):
    """
    figures out what the color is of the cell
    :param code: location code of cell
    :param gen: generation of cell used to calculate size of cell
    :param grid: grid that "contains" cell
    :param problem_size: max amount of cell divisions
    :return: white = "W", gray = "G", black = "B"
    """
    cell_x = int('{0:06b}'.format(code)[1::2], 2)  # extract interlaced x coordinate
    cell_y = int('{0:06b}'.format(code)[0::2], 2)  # extract interlaced y coordinate
    # print 'checking {0:03b} , {1:03b} from {2:06b} in gen {3:d}'.format(cell_x, cell_y, code, gen)
    score = 0
    cell_length = (2 ** (problem_size - gen))
    # print 'cell size is {0:d} x {0:d}'.format(cell_length)
    for i in range(cell_length):
        for j in range(cell_length):
            score += grid[cell_y + i][cell_x + j]

    # print 'score of {0:s} is {1:d}'.format(binary_to_quaternary_string(code, problem_size), score)
    if score == 0:
        return 'W'
    elif score == cell_length ** 2:  # 2 ** (problem_size - gen) means entire cell was black
        return 'B'
    else:
        return 'G'


def tree_as_string(tree, problem_size):
    """
    dump contents of quadtree in prettified format
    :param tree:
    :param problem_size:
    :return:
    """
    desired = '(\n'
    for cell in tree:
        cell = list(cell)
        for i in range(len(cell)):
            if isinstance(cell[i], int):
                if i == 0:
                    cell[i] = binary_to_quaternary_string(cell[i], problem_size)
                elif cell[i] == sys.maxsize:
                    cell[i] = '#'
                elif cell[i] == sys.maxsize + 1:
                    cell[i] == 'ERROR'
                else:
                    cell[i] = '%d' % cell[i]
        desired += str(cell) + ',\n'
    desired += ')'
    return desired


def divide(cell, grid, problem_size):
    for i in range(4):
        if cell[3 + i] != sys.maxsize:
            cell[3 + i] -= 1
    generations = [cell[1] + 1] * 4
    new_code_bits = map(lambda x: x << (2 * (problem_size - generations[x])), range(4))
    codes = map(lambda x: cell[0] | x, new_code_bits)
    # print 'Children codes are: %s' % str(map(lambda x: binary_to_quaternary_string(x, problem_size), codes))
    colors = map(lambda x: fetch_color(codes[x], generations[x], grid, problem_size), range(4))
    east_levels = [0, cell[3], 0, cell[3]]
    north_levels = [0, 0, cell[4], cell[4]]
    west_levels = [cell[5], 0, cell[5], 0]
    south_levels = [cell[6], cell[6], 0, 0]
    return map(lambda x: list(x), zip(codes, generations, colors, east_levels, north_levels, west_levels, south_levels))


def populate_tree(grid, lqtld):
    """
    populate the linear quad tree with cells that track level differences.
    :param grid: input data that needs to be represented as quadtree
    :param lqtld: list structure that will be filled with cells
    :param current_gen: current iteration of cell divisions
    :return: this is a non-returning function
    """
    problem_size = int(math.log(len(grid), 2))  # determine max amount of cell divisions

    while grey_nodes_in_tree(lqtld):  # line 3 of figure 7 from paper
        # print 'current step: %s' % tree_as_string(copy.deepcopy(lqtld), problem_size)
        for i in range(len(lqtld)):
            if lqtld[i][2] == 'G':  # find first gray node
                # print 'found gray node with code: {0:s} and gen: {1:d}'.format(binary_to_quaternary_string(lqtld[
                # i][0], problem_size), lqtld[i][1])
                update_neighbors(lqtld, lqtld[i], problem_size)  # line 5 of figure 7
                break
        target_node = lqtld.pop(i)
        children = divide(target_node, grid, problem_size)
        for child in children:
            update_neighbors(lqtld, child, problem_size)
        lqtld.extend(children)


if __name__ == '__main__':
    flipped_grid = [[1, 1, 1, 1, 1, 0, 0, 0],  # this is sample data from paper
                    [1, 1, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
    desired_grid = np.flipud(np.mat(flipped_grid))
    linear_tree = [[0, 0, 'G', sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize]]
    populate_tree(desired_grid.tolist(), linear_tree)
    linear_tree = sorted(linear_tree, key=lambda x: x[1], reverse=False)
    # print tree_as_string(linear_tree, int(math.log(len(desired_grid.tolist()), 2)))
