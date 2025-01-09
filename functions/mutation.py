"""
The module contains a function 'mutation' that performs mutation operation
"""

from ase import Atoms, Atom
import math
import random
import numpy as np
from functions.small_functions import get_r, calc_distance, get_rand_xyz, get_image_positions

def mutation(atoms):
    """
    Performs mutation of the given structure by changing the position of one atom between the layers.

    Args:
        atoms (ase.Atoms): Structure to be mutated.

    Returns:
        ase.Atoms: The structure created by mutation.
    """
    # An output structure that is a copy of the one given as an argument to the function
    structure = atoms.copy()
    # Symbol of the atom located between the layers
    atom_symbol = structure[len(structure)-1].symbol

    # Dimensions of the unit cell
    cell = structure.get_cell()
    array = structure.cell.cellpar()
    a = array[0] # Dimension in the x-axis direction
    b = array[1]
    b_y = b * np.sqrt(3) / 2  # Component b in the y direction
    b_x = b / 2  # Component b in the x direction
    c = array[2] # Dimension in the z-axis direction
    c_half = c / 2 # Half the height of the cell

    atom_indexes = [] #  List with indexes of atoms lying between layers

    # Adding indexes of atoms lying between layers to the list of atom_indexes
    for atom in structure:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            atom_indexes.append(atom.index)

    index_to_del = random.choice(atom_indexes) # Randomly selected atom index to be removed
    del structure[index_to_del] # Removing an atom from the structure
    atom_indexes.remove(index_to_del) # Removing the index of the removed atom from the list

    # Temporary structure storing the positions of atoms between layers
    # needed to create the image_structures object
    tmp_structure = Atoms()
    tmp_structure.set_cell(cell)
    for atom in structure:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            tmp_structure.append(atom)

    tol_r = 0.1 # Tolerance used when checking the distance between atoms

    # Temporary 3x3 structure that stores the positions of atoms between the layers and their images
    image_structure = Atoms()
    image_structure.set_cell(cell)
    image_structure = image_structure * [3, 3, 1]

    # Adding atoms to image positions based on the tmp_structures object
    for atom in tmp_structure:
        image_positions = get_image_positions(cell, atom)

        for position in image_positions:
            image_structure.append(Atom(atom_symbol, position))

    # Flag informing whether it was possible to select an atom position that does not collide
    # with the atomic images and atoms from layers in the target structure
    should_continue = True

    while should_continue:
        collision_with_image = False # Flag informing about collisions with atomic images
        random_positions = get_rand_xyz(cell)
        random_atom = Atom(atom_symbol, random_positions) # Atom with new, random coordinates
        middle_positions = random_positions + np.array([a - b_x, b_y, 0.])
        middle_atom = Atom(atom_symbol, middle_positions) # Atom with random coordinates moved to the center of image_structure

        for atom in image_structure: # Checking the distance between middle_atom and images
            d = calc_distance(middle_atom, atom)
            if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_image = True
                    break

        if not collision_with_image:
            collision_with_structure = False # Flag informing about a collision with atoms in the target structure
            for atom in structure: # Checking the distance between random_atom and atoms in the target structure
                d = calc_distance(random_atom, atom)
                if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_structure = True
                    break

            # If there was no collision - add an atom to the target structure
            if not collision_with_structure:
                structure.append(random_atom)
                should_continue = False

    return structure
