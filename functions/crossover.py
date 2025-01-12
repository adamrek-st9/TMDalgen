"""
The module contains a function 'crossover' that performs crossover operation.
"""

from ase import Atoms
import math
import random
from itertools import combinations
from functions.small_functions import get_r, second_part
from functions.prep_struct import prep_struct

def crossover(atoms1, atoms2, n, struct_filename, struct_size):
    """
    Performs crossover between two given structures by exchanging atoms between them.

    Args:
        atoms1 (ase.Atoms): First structure.
        atoms2 (ase.Atoms): Second structure
        n (int): Number of atoms exchanged between structures during crossover.
        struct_filename (str): Filename of basic dichalcogenide structure (e.g., MoS2).
        struct_size (str): Size of structures (e.g., 4x4).

    Returns:
        ase.Atoms: The structure created by crossover.
    """
    # Parent structures - the same cells as atoms1 and atoms2
    # Atoms located between the layers in atoms1 and atoms2 will be added to them
    parent1 = Atoms()
    parent1.set_cell(atoms1.cell)
    parent2 = Atoms()
    parent2.set_cell(atoms2.cell)

    # Output child object
    child = prep_struct(struct_filename, struct_size)

    #  Dimensions of the parents unit cell
    cell = atoms1.get_cell()
    array = atoms1.cell.cellpar()
    c = array[2]  # Dimension in the direction of the z-axis
    c_half = c / 2 # Half the height of the cell

    tol_r = 0.1  # Tolerance used when checking the distance between atoms

    # Adding atoms located between layers in atoms1 and atoms2 structures to the parent structures
    for atom in atoms1:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            parent1.append(atom)
    for atom in atoms2:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            parent2.append(atom)

    combinations1 = list(combinations(parent1, n)) # All n-atomic combinations from parent1
    combinations2 = list(combinations(parent2, n)) # All n-atomic combinations from parent2

    while True:
        # Temporary structures of children
        tmp_child1 = Atoms()
        tmp_child1.set_cell(cell)
        tmp_child2 = Atoms()
        tmp_child2.set_cell(cell)
        # Part 'a' of parent1 made of atoms selected during combination
        parent1_a = Atoms(random.choice(combinations1))
        # Part 'b' of parent1 made of remaining atoms
        parent1_b = second_part(parent1, parent1_a)
        # Part 'a' of parent2 made of atoms selected during combination
        parent2_a = Atoms(random.choice(combinations2))
        # Part 'b' of parent2 made of remaining atoms
        parent2_b = second_part(parent2, parent2_a)

        # Adding atoms between layers to child no. 1
        for atom in parent1_a + parent2_b:
            tmp_child1.append(atom)

        # Adding atoms between layers to child no. 2
        for atom in parent1_b + parent2_a:
            tmp_child2.append(atom)

        # Collision check between atoms and their images in tmp_child1
        image_structure1 = tmp_child1 * [3, 3, 1]
        collision_in_child1 = False
        for i in range(len(image_structure1)):
            for j in range(i + 1, len(image_structure1)):
                distance = image_structure1.get_distance(i, j)
                if distance < (2 * get_r(image_structure1.get_atomic_numbers()[0]) + tol_r):
                    collision_in_child1 = True
                    break

        # Collision check between atoms and their images in tmp_child2
        image_structure2= tmp_child2 * [3, 3, 1]
        collision_in_child2 = False
        for i in range(len(image_structure2)):
            for j in range(i + 1, len(image_structure2)):
                distance = image_structure2.get_distance(i, j)
                if distance < (2 * get_r(image_structure2.get_atomic_numbers()[0]) + tol_r):
                    collision_in_child2 = True
                    break

        # If there is no collision in tmp_child1 or tmp_child2 - create output child structure
        if not collision_in_child1:
            for atom in tmp_child1:
                child.append(atom)
            break
        if not collision_in_child2:
            for atom in tmp_child2:
                child.append(atom)
            break

    return child
