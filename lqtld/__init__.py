import sys
import math
import numpy as np


# Node format: (code number, generation, color, east, north, west, south)
# NOTE: sys.maxsize used to denote presence of wall as neighbor
# NOTE: 0 = bottom left, 1 = bottom right, 2 = top left, 3 = top right
problem_size = 3
def qlao(code, tx, ty, direction):
    """
    This is the circled plus operator from the paper. QLAO stands for quad location addition operator
    :rtype: int
    :param code:
    :param tx: the tx value fom the paper
    :param ty: the ty value from the paper
    :param direction: direction vector in interleaved binary format
    :return:
    """
    return (((code | ty) + (direction & tx)) & tx) | (((code | tx) + (direction & ty)) & ty)


def update_neighbors(tree, cell):
    """
    updates of all neighbors of a given cell code regardless of brotherhood. the first for loop of figure 7
    :param tree: linear tree containing neighbors
    :param cell: tuple representing gray cell about to be divided
    :return: this is a non-returning functions
    """
    # generate direction increments for finding neighbors
    # NOTE: directions in list are in same order as in tuple, this is used in the loop
    directions = [int('01', 2), int('10', 2), int('01' * problem_size, 2),
                  int('10' * problem_size, 2)]  # west, south, north, east
    neighbor_codes = []
    directions = list(map(lambda x: x << (2 * (problem_size - cell[1])), directions))
    for i in range(4):
        if cell[3 + i] == sys.maxsize:  # if neighbor is known wall, make no changes to any LD
            # print 'Neighbor {0:d} is wall, doing nothing'.format(i)
            neighbor_codes.append(sys.maxsize)
        else:
            neighbor_codes.append(qlao(cell[0], int('01' * problem_size, 2), int('10' * problem_size, 2),
                                       directions[i]))  # this is theorem 1 from paper
            neighbor = next((x for x in tree if x[0] == neighbor_codes[i]), None)
            if neighbor is not None and neighbor[1] == cell[1]:
                neighbor[3 + ((2 + i) % 4)] += 1


def grey_nodes_in_tree(tree):
    """
    check if there are gray nodes in the tree using built-in next function
    :param tree:
    :return: whether or not there's a gray node in the tree
    """
    return next((node for node in tree if node[2] == 'G'), None) is not None


def binary_to_quaternary_string(code):
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


def fetch_color(code, gen, grid):
    """
    figures out what the color is of the cell
    :param code: location code of cell
    :param gen: generation of cell used to calculate size of cell
    :param grid: grid that "contains" cell
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


def tree_as_string(tree):
    """
    dump contents of quadtree in prettified format
    :param tree:
    :return:
    """
    desired = '(\n'
    for cell in tree:
        cell = list(cell)
        for i in range(len(cell)):
            if isinstance(cell[i], int):
                if i == 0:
                    cell[i] = binary_to_quaternary_string(cell[i])
                elif cell[i] == sys.maxsize:
                    cell[i] = '#'
                elif cell[i] == sys.maxsize + 1:
                    cell[i] == 'ERROR'
                else:
                    cell[i] = '%d' % cell[i]
        desired += str(cell) + ',\n'
    desired += ')'
    return desired


def divide(cell, grid):
    for i in range(4):
        if cell[3 + i] != sys.maxsize:
            cell[3 + i] -= 1
    generations = [cell[1] + 1] * 4
    new_code_bits = map(lambda x: x << (2 * (problem_size - generations[x])), range(4))
    codes = list(map(lambda x: cell[0] | x, new_code_bits))
    # print 'Children codes are: %s' % str(map(lambda x: binary_to_quaternary_string(x, problem_size), codes))
    colors = list(map(lambda x: fetch_color(codes[x], generations[x], grid), range(4)))
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
                update_neighbors(lqtld, lqtld[i])  # line 5 of figure 7
                break
        target_node = lqtld.pop(i)  # pop gray node off
        children = divide(target_node, grid)
        for child in children:
            update_neighbors(lqtld, child)  # check non-brother neighbors of kids
        lqtld.extend(children)  # add new kids

def get_neighbor(cell, level_dif_index, direction):
    """
    returns neighbor code in given direction
    :param cell: cell whose neighbor we want
    :param level_dif_index: level difference of neighbor in that direction
    :param direction: direction vector to apply to the code
    :param problem_size: maximum amount of cell divisions
    :return: location code of neighbor
    """
    dd = cell[level_dif_index]
    if dd != sys.maxsize:  # if the direction requested is not facing wall
        n_q, l = cell[:2]
        shift = 2 * (problem_size - l - dd)
        tx = int('01' * problem_size, 2)
        ty = int('10' * problem_size, 2)
        if dd < 0:
            return qlao((n_q >> shift << shift), tx, ty, (direction << shift))
        else:
            return qlao(cell[0], tx, ty, direction << (2 * (problem_size - l)))

def find_containing_cell(row, col, tree):
    row_bits = '{0:03b}'.format(row)
    col_bits = '{0:03b}'.format(col)
    pixel_code = int(''.join([row + col for row, col in zip(row_bits, col_bits)]), 2)
    query_code = -1
    print 'Pixel Code is:' + binary_to_quaternary_string(pixel_code, 3) + ' (' + col_bits + ', ' + row_bits + ')'
    for i in range(1, problem_size):
        mask = int('1'*(2*i) + '0'*(2 * (problem_size - i)), 2)
        print bin(pixel_code) + ' & ' + bin(mask)
        query_code = pixel_code & mask
        query_gen = i
        print 'Searching for cell with code: {0:06b}'.format(query_code) + ' and gen: ' + str(query_gen)
        result = binary_search(tree, query_code, query_gen)
        if result != -1:
            break

    return 'Containing Cell: ' + binary_to_quaternary_string(query_code, 3)



def linear_search(code, gen, tree):
    print 'Searching for code: ' + binary_to_quaternary_string(code) + ' of gen: ' + str(gen)
    for i in range(len(tree)):
        if code == tree[i][0] and gen == tree[i][1]:
            return i
    return -1

def get_all_neighbor_indices(index, tree):
    print 'getting neighbor_indices of: ' + binary_to_quaternary_string(tree[index][0])
    direction_vectors = [1, 2, int('01'*problem_size,2), int('10'*problem_size, 2)]  # generate diretion vectors
    neighbor_indices = []
    for i in range(len(direction_vectors)):
        if tree[index][3+i] != sys.maxsize: #if there isn't a wall in that direction
            neighbor_code = get_neighbor(tree[index], 3+i, direction_vectors[i])
            neighbor_gen = tree[index][1]+tree[index][3+i]
            print 'Checking neighbor: ' + binary_to_quaternary_string(neighbor_code) + ' of gen: ' + str(neighbor_gen)
            neighbor_index = linear_search(neighbor_code, neighbor_gen, tree)
            if tree[index][3+i] <=0: #neighbor is bigger or same size, just add it on
                neighbor_indices.append(linear_search(neighbor_code, neighbor_gen, tree))
            else:  # otherwise, you're gonna have to find all the cells flush to that side
                get_sub_cells(neighbor_code, tree[index][1], i, neighbor_indices, tree)
    return neighbor_indices

def get_sub_cells(code, gen, direction, index_list, tree):
    """
    load all subcells flush to a certain side into a list
    """
    if gen == problem_size:
        print 'Couldn\'t find cell' + binary_to_quaternary_string(code) + ' of gen: ' + str(gen)
    flush_combos = [[0,2],[0,1],[1,3],[2,3]] #new bits to add given each case
    codes_to_check = map(lambda x: code | (x << (2 * (problem_size - gen - 1))), flush_combos[direction])  # create child codes
    for code in codes_to_check:
        child_index = linear_search(code, gen+1, tree)
        if child_index != -1:
            index_list.append(child_index)
        else:
            get_sub_cells(code, gen+1, direction, index_list, tree)

if __name__ == '__main__':
    flipped_grid = [[1, 1, 1, 1, 1, 0, 0, 0],  # this is sample data from paper
                    [1, 1, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
    desired_grid = np.flipud(np.mat(flipped_grid))  # flip my matrix to fit coordinates of code system
    linear_tree = [[0, 0, 'G', sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize]]  # instantiate tree
    populate_tree(desired_grid.tolist(), linear_tree)  # do the thing with the stuff
    linear_tree = sorted(linear_tree, key=lambda x: x[0], reverse=False)  # sort tree to match paper
    print tree_as_string(linear_tree)
    test_neighbors = []
    for index in get_all_neighbor_indices(5, linear_tree):
        test_neighbors.append(linear_tree[index])
    print tree_as_string(test_neighbors)
