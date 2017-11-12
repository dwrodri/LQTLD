import math
import sys
import png
import numpy as np
class LQTLD:
    def __init__(self, occupancy_matrix, map_resolution, dimensions, filename=None):
        self.occupancy_matrix = occupancy_matrix
        self.map_resolution = map_resolution
        self.width, self.height = dimensions
        if self.width != self.height or not filter(lambda x: x & (x-1) == 0, dimensions):  # if not a power of two and not square
             self.pad(int(math.log(len(occupancy_matrix), 2))+1)
             self.width = self.height = len(occupancy_matrix)
        self.r = int(math.log(len(occupancy_matrix), 2))
        self.tree = []
        if filename is None:
            self.tree = [[0, 0, 'G', sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize]]  # initialize tree with node representing entire map
            self.generate_tree(self.tree)  # this is where the magic happens...
            self.tree = sorted(self.tree, key=lambda x: x[0])
            self.write_tree_to_file()
        else:
            with open(filename) as tree_data:
                for line in tree_data:
                    cell = []
                    for element in line.split() :
                        if element != 'W' and element != 'B':
                            cell.append(int(element))
                        else:
                            cell.append(element)

                    self.tree.append(cell)
        print self.r
        print self.repr_tree()

    def write_tree_to_file(self):
        file_object = open('/tmp/test_lqtld.txt', 'w')
        for cell in self.tree:
            line = ''
            for value in cell:
                line += str(value) + ' '
            file_object.write(line + '\n')

    def pad(self, new_r):
        """
        add obstacle values into the occupancy_matrix in order to make both dimensions the same power of two
        :param new_r: the new order of magnitude for the matrix
        """
        horiz_extension = (2**new_r)-len(self.occupancy_matrix[0])
        vert_extension = (2**new_r)-len(self.occupancy_matrix)
        for row in self.occupancy_matrix:  # pad existing rows
            row.extend([1]*horiz_extension)
        for i in range(vert_extension):
            self.occupancy_matrix.append([1]*(2**new_r))

    def generate_tree(self, grid):
        """
        populate the linear quad tree with cells that track level differences.
        :param self.tree: list structure that will be filled with cells
        :param current_gen: current iteration of cell divisions
        :return: this is a non-returning function
        """
        while self.grey_nodes_in_tree():  # lqtld includes GRAY areas do
            for i in range(len(self.tree)):
                if self.tree[i][2] == 'G':  # find first gray node
                    self.update_neighbors(self.tree[i])  # neighbors of first GRAY node increased by one
                    break
            target_node = self.tree.pop(i)  # delete it
            children = self.divide(target_node)  # append its four children
            for child in children:  # for each equal-sized neighbor of each child (that is not brother)
                self.update_neighbors(child)  # corresponding level differences increased by 1
            self.tree.extend(children)  # add new kids since they were out of the list

    def update_neighbors(self, cell):
        """
        updates of all neighbors of a given cell code regardless of brotherhood. the first for loop of figure 7
        :param tree: linear tree containing neighbors
        :param cell: tuple representing gray cell about to be divided
        :param self.r: problem size used to generate direction codes
        :return: this is a non-returning functions
        """
        # generate direction increments for finding neighbors
        # NOTE: directions in list are in same order as in tuple, this is used in the loop
        directions = [int('01', 2), int('10', 2), int('01' * self.r, 2),
                      int('10' * self.r, 2)]  # east, north, west, south
        neighbor_codes = []
        directions = map(lambda x: x << (2 * (self.r - cell[1])), directions)
        for i in range(4):
            if cell[3 + i] == sys.maxsize:  # if neighbor is known wall, make no changes to any LD value
                neighbor_codes.append(sys.maxsize)
            else:
                neighbor_codes.append(self.qlao(cell[0], int('01' * self.r, 2), int('10' * self.r, 2), directions[i]))  # this is theorem 1 from paper
                neighbor = next((x for x in self.tree if x[0] == neighbor_codes[i]), None)
                if neighbor is not None and neighbor[1] == cell[1]:
                    neighbor[3 + ((2 + i) % 4)] += 1

    def divide(self, cell):
        """
        splits a cell into its four children
        :param cell: the cell to be divided
        :return: list of four children
        """
        for i in range(4):
            if cell[3 + i] != sys.maxsize:
                cell[3 + i] -= 1
        generations = [cell[1] + 1] * 4
        new_code_bits = map(lambda x: x << (2 * (self.r - generations[x])), range(4))
        codes = map(lambda x: cell[0] | x, new_code_bits)
        colors = map(lambda x: self.fetch_color(codes[x], generations[x]), range(4))  # fetch colors for all new codes
        east_levels = [0, cell[3], 0, cell[3]]  # east levels go to SE and NE corners
        north_levels = [0, 0, cell[4], cell[4]]  # north levels go NE and NW corners
        west_levels = [cell[5], 0, cell[5], 0]  # west levels go to SW and SW corners
        south_levels = [cell[6], cell[6], 0, 0]  # south levels go to SW and SE corners
        return map(lambda x: list(x), zip(codes, generations, colors, east_levels, north_levels, west_levels, south_levels))  # zip and convert kids to lists

    def code_to_pixel(self, code):
        """
        convert cell address code to pixels.
        :param code: the cell code to be converted
        :return: list of [row, col] that can be used by other programs
        """
        y_bits = bin(code).zfill(2*self.r)[::2]  # print length of binary matching code
        x_bits = bin(code).zfill(2*self.r)[1::2]
        row = 0
        col = 0
        for i in range(self.r):
            x_bit = x_bits[i]  # get ith bit from each one
            y_bit = y_bits[i]
            if x_bit == '1':  # if bit is one...
                col += len(self.occupancy_matrix[0]) >> (i+1)  # add length of cell generation
            if y_bit == '1':
                row += len(self.occupancy_matrix) >> (i+1)  # do same for y bits
        #print 'locaiton on matrix is: %d, %d' % (row, col)
        return [row, col]


    def pose_to_pixel(self, pose_coords):
        """
        converts x,y on ros stage to row, col on occupancy_matrix
        :param pose_coords: x,y in iterable format
        """
        pose_x, pose_y = pose_coords  # unpack coordinates
        pixel_row = int((pose_y/self.map_resolution) + (len(self.occupancy_matrix)-224)/2)
        pixel_col = int((pose_x/self.map_resolution) + (len(self.occupancy_matrix[0])-24)/2)
        return [pixel_row, pixel_col]

    def pixel_to_pose(self, pixel_coords):
        """
        converts row, col to ros stage x, y
        :param pixel_coords: row, col in iterable format
        """
        pixel_row, pixel_col = pixle_coords
        pose_x = (pixel_col - (len(self.occupancy_matrix[0])-24)/2) * self.map_resolution
        pose_y = (pixel_row - (len(self.occupancy_matrix)-224)/2) * self.map_resolution
        return [pose_x, pose_y]

    def get_all_neighbor_indices(self, index):
        direction_vectors = [1, 2, int('01'*self.r,2), int('10'*self.r, 2)]  # generate diretion vectors
        neighbor_indices = []
        for i in range(len(direction_vectors)):
            if self.tree[index][3+i] != sys.maxsize: #if there isn't a wall in that direction
                neighbor_code = self.get_neighbor(self.tree[index], 3+i, direction_vectors[i])
                neighbor_gen = self.tree[index][1]+self.tree[index][3+i]
                neighbor_index = self.linear_search(neighbor_code, neighbor_gen)
                if self.tree[index][3+i] <=0: #neighbor is bigger or same size, just add it on
                    neighbor_indices.append(self.linear_search(neighbor_code, neighbor_gen))
                else:  # otherwise, you're gonna have to find all the cells flush to that side
                    self.get_sub_cells(neighbor_code, self.tree[index][1], i, neighbor_indices)
        return neighbor_indices

    def get_sub_cells(self, code, gen, direction, index_list):
        """
        load all subcells flush to a certain side into a list
        """
        if gen == self.r:
            print 'Couldn\'t find cell' + binary_to_quaternary_string(code) + ' of gen: ' + str(gen)
        flush_combos = [[0,2],[0,1],[1,3],[2,3]] #new bits to add given each case
        codes_to_check = map(lambda x: code | (x << (2 * (self.r - gen - 1))), flush_combos[direction])  # create child codes
        for code in codes_to_check:
            child_index = self.linear_search(code, gen+1)
            if child_index != -1:
                index_list.append(child_index)
            else:
                self.get_sub_cells(code, gen+1, direction, index_list)

    def get_open_cells(self):
        """
        return a list of only the white cells
        """
        return filter(lambda x: x[2] == 'W', self.tree)

    def fetch_color(self, code, gen):
        """
        figures out what the color is of the cell
        :param code: location code of cell
        :param gen: generation of cell used to calculate size of cell
        :param grid: grid that "contains" cell
        :param self.r: max amount of cell divisions
        :return: white = "W", gray = "G", black = "B"
        """
        cell_row, cell_col = self.code_to_pixel(code)  # get row and col of cell
        cell_size = 2 ** (self.r - gen)
        score = 0
        #print 'The cell being evaluated is code' + self.binary_to_quaternary_string(code)
        for i in range(cell_size):
            for j in range(cell_size):
                score += self.occupancy_matrix[cell_row + i][cell_col + j]  # black cells are represented with 1s on the grid

        if score == 0:
            return 'W'
        elif score == cell_size:
            return 'B'
        else:
            return 'G'

    def repr_tree(self):
        """
        pretty print for tree
        """
        tree_string = '[\n'
        for cell in self.tree:
            cell[0] = self.binary_to_quaternary_string(cell[0])
            tree_string += str(cell) + ',\n'
            cell[0] = int(cell[0], 4)
        tree_string += ']'
        return tree_string

    def binary_to_quaternary_string(self, code):
        """
        convert the int location code to printable quaternary value
        :param code: binary form of code to be printed
        :param self.r: amount of quaternary values
        :return: string form of quaternary
        """
        if code == sys.maxsize:
            return '#'
        desired = ''
        code_bits = bin(code)[2:].zfill(2*self.r)
        for i in range(1, 2*self.r+1, 2):
            desired += str(int(code_bits[i-1:i+1], 2))
        return desired

    def grey_nodes_in_tree(self):
        """
        check if there are gray nodes in the tree using built-in next function
        :param tree:
        :return: whether or not there's a gray node in the tree
        """
        return next((node for node in self.tree if node[2] == 'G'), None) is not None

    def qlao(self, code, tx, ty, direction):
        """
        This is the circled plus operator from the paper. QLAO stands for quad location addition operator
        :rtype: int
        :param code:
        :param tx: the tx value fom the paper
        :param ty: the ty value from the paper
        :param direction: direction vector in interleaved binary format
        :return: the calculated neighborv code based on this operation
        """
        return (((code | ty) + (direction & tx)) & tx) | (((code | tx) + (direction & ty)) & ty)

    def get_neighbor(self, cell, level_dif_index, direction):
        """
        returns neighbor code in given direction
        :param cell: cell whose neighbor we want
        :param level_dif_index: level difference of neighbor in that direction
        :param direction: direction vector to apply to the code
        :return: location code of neighbor
        """
        dd = cell[level_dif_index]  # this is the level differential
        if dd != sys.maxsize:  # if the direction requested is not facing wall
            n_q, l = cell[:2]  # n_q is the neighbor code, l is level
            shift = 2 * ( self.r - l - dd)  # this is the amuont of bitshifting required to comqnesate for level differences
            tx = int('01' *  self.r, 2)  # bitmasks used in QLAO
            ty = int('10' *  self.r, 2)
            if dd < 0:  # if neighbor is higher up in tree
                return self.qlao((n_q >> shift << shift), tx, ty, (direction << shift))
            else:  # if neighbor is further down or at same level...
                return self.qlao(cell[0], tx, ty, direction << (2 * ( self.r - l)))

    def shit_get(self, row, col):
        """
        pls bb no
        """
        row_bits = bin(row)[2:].zfill(self.r)
        col_bits = bin(col)[2:].zfill(self.r)
        pixel_code = int(''.join([row + col for row, col in zip(row_bits, col_bits)]), 2)
        query_code = -1
        result = -1
        print 'Pixel Code is:' + self.binary_to_quaternary_string(pixel_code) + ' (' + col_bits + ', ' + row_bits + ')'
        for i in range(1, self.r):
            mask = int('1'*(2*i) + '0'*(2 * (self.r - i)), 2)
            print bin(pixel_code) + ' & ' + bin(mask)
            query_code = pixel_code & mask
            query_gen = i
            print 'Searching for cell with code:' + self.binary_to_quaternary_string(query_code) + ' and gen: ' + str(query_gen)
            result = self.linear_search(query_code, query_gen)
            if result != -1:
                break

        return result

    def linear_search(self, code, gen):
        for i in range(len(self.tree)):
            if code == self.tree[i][0] and gen == self.tree[i][1]:
                return i
        return -1

    def get_center_pixel_from_index(self, index):
        edge_row, edge_col = self.code_to_pixel(self.tree[index][0])
        half = 2**(self.r - self.tree[index][1] - 1)
        return [edge_row+half, edge_col+half]


    def color_cell(self, cell_index, color):
        cell = self.tree[cell_index]
        cell_row, cell_col = self.code_to_pixel(cell[0])  # get corner of cell
        for i in range(len(self.occupancy_matrix)>>cell[1]):
            for j in range(len(self.occupancy_matrix[0])>>cell[1]):
                if i == 0 or j == 0 or i == (len(self.occupancy_matrix)>>cell[1])-1 or j == (len(self.occupancy_matrix[0])>>cell[1])-1:
                    self.occupancy_matrix[cell_row+i][cell_col+j] = color

    def draw_all_usable_cells(self):
        """
        edits the occupancy_matrix to draw in the edges of the cells
        """
        for i in range(len(self.tree)):
            if self.tree[i][2] != 'B' and self.tree[i][1] < (self.r-1):
                self.color_cell(i, int('0x8F8F8F', 16))

    def draw_pretty_quarters(self):
        for i in range(len(self.tree)):
            if self.tree[i][1] < (self.r-1):
                header = int(self.tree[i][0]>>18)
                if header == 0:
                    self.color_cell(i, int('0xFF0000', 16))  # 0 is red
                elif header == 1:
                    self.color_cell(i, int('0x008000', 16))  # 1 is lime
                elif header == 2:
                    self.color_cell(i, int('0x0000FF', 16))  # 2 is Blue
                else:
                    self.color_cell(i, int('0x800080', 16))  # 3 is dark purple

    def draw_point(self, point_row, point_col, color, size):
        for i in range(size):
            for j in range(size):
                self.occupancy_matrix[point_row+i][point_col+j] = color
                self.occupancy_matrix[point_row+i][point_col-j] = color
                self.occupancy_matrix[point_row-i][point_col+j] = color
                self.occupancy_matrix[point_row-i][point_col-j] = color

    def generate_debug_png(self, filename):
        """
        debugger function that uses PyPNG to render an RGB image used to debug stuff
        :param filename: string used as filename for output_PNG
        """
        output_file = open(filename, 'wb')
        writer = png.Writer(2**self.r, 2**self.r)
        pixel_matrix  = []
        for row in range(len(self.occupancy_matrix)):
            pixel_matrix.append([])
            for col in range(len(self.occupancy_matrix[0])):
                if self.occupancy_matrix[row][col] == 0:
                    pixel_matrix[row].extend([255,255,255])
                elif self.occupancy_matrix[row][col] == 1:
                    pixel_matrix[row].extend([0,0,0])
                else:
                    b = self.occupancy_matrix[row][col] & 255
                    g = (self.occupancy_matrix[row][col]>>8) & 255
                    r = (self.occupancy_matrix[row][col]>>16) & 255
                    pixel_matrix[row].extend([r,g,b])
        writer.write(output_file, pixel_matrix[::-1])

if __name__ == '__main__':
    flipped_grid = [[1, 1, 1, 1, 1, 0, 0, 0],  # this is sample data from paper
                    [1, 1, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
    desired_grid = np.flipud(np.mat(flipped_grid)).tolist()  # flip my matrix to fit coordinates of code system
    lqtld = LQTLD(desired_grid, 1, [8, 8])
    lqtld.draw_all_usable_cells()
    lqtld.generate_debug_png('/tmp/mini_test.png')
