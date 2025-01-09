"""
The module contains a function 'gen_rand_struct' that generates a random structure.
"""

from ase import Atoms, Atom
import numpy as np
from functions.small_functions import get_r, calc_distance, get_image_positions, get_rand_xyz
from functions.prep_struct import prep_struct

def gen_rand_struct(structure_file_name, size, atom_symbol, n_atoms):
    """
    Generates a two-layer structure with atoms randomly distributed between the layers, based on the given parameters.

    Args:
        structure_file_name (str): Name of the file containing the dichalcogenide structure.
        size (str): Size of the structure (e.g., 4x4).
        atom_symbol (str): Chemical symbol of atoms between the layers.
        n_atoms (int): Number of atoms between layers.

    Returns:
        ase.Atoms: The generated structure.
    """
    # The output structure to which atoms will be added
    structure = prep_struct(structure_file_name, size)

    # Dimensions of the unit cell
    cell = structure.get_cell()
    array = cell.cellpar()
    a = array[0] # Dimension in the x-axis direction
    b = array[1]
    b_y = b * np.sqrt(3) / 2 # Component b in the y direction
    b_x = b / 2 # Component b in the x direction

    # Temporary structure storing the positions of atoms between layers
    tmp_structure = Atoms()
    tmp_structure.set_cell(cell)

    # Temporary 3x3 structure that stores the positions of atoms between the layers and their imagess
    image_structure = Atoms()
    image_structure.set_cell(cell)
    image_structure = image_structure * [3, 3, 1]

    tol_r = 0.1 # Tolerance used when checking the distance between atoms

    # Adding atoms in a given number
    for i in range(n_atoms):
        while True:
            # Atom with random coordinates
            random_positions = get_rand_xyz(cell)
            random_atom = Atom(atom_symbol, random_positions)

            # Atom with random coordinates shifted by vector a and b
            # Middle atom in a 3x3 temporary structure
            middle_positions = random_positions + np.array([a - b_x, b_y, 0.])
            middle_atom = Atom(atom_symbol, middle_positions)

            # Positions to which atoms will be added in the temporary 3x3 structure
            image_positions = get_image_positions(cell, random_atom)

            collision_with_atom = False # Flag informing about collision with atoms from layers
            # Checking the distance between the added atom and the atoms from the layers
            for atom in structure:
                d = calc_distance(random_atom, atom)
                if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_atom = True
                    break

            if not collision_with_atom:
                collision_with_image_atom = False # Flag informing about collisions with atomic images
                # Checking whether the added atom does not collide with the images
                for atom in image_structure:
                    d = calc_distance(middle_atom, atom)
                    if d < (get_r(middle_atom.number) + get_r(atom.number) + tol_r):
                        collision_with_image_atom = True
                        break

                if not collision_with_image_atom:
                    # Adding atoms to the structure with images in all possible positions
                    for position in image_positions:
                        image_structure.append(Atom(atom_symbol, position))
                    # Adding an atom that does not interfere with layers and images to the target structure
                    structure.append(random_atom)
                    break

    return structure
