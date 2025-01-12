"""
The module contains a function 'get_r', 'calc_distance', 'second_part', 'get_rand_xyz' and 'get_image_positions' -
small functions used in other modules.
"""

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
import random

def get_r(atomic_number):
    """
    Calculate the radius of an atom.

    Args:
        atomic_number (int): The atomic number.

    Returns:
        float: The radius of an atom.
    """
    tmp_r = 0.9 * covalent_radii[atomic_number]
    return tmp_r

def calc_distance(atom1, atom2):
    """
    Calculate the distance between two atoms.

    Args:
        atom1 (ase.Atom): The first atom.
        atom2 (ase.Atom): The second atom.

    Returns:
        float: The distance between atom1 and atom2.
    """
    tmp_d = np.linalg.norm(atom1.position - atom2.position)
    return tmp_d

def second_part(parent, part):
    """
    Returns atoms from the parent structure that are not in the part object.
    Function used in the crossover module.

    Args:
        parent (ase.Atoms): The parent structure.
        part (ase.Atoms): The part object.

    Returns:
        ase.Atoms: The structure with atoms from the parent object that are not in the part object.
    """
    tmp_second_part = Atoms()
    tmp_second_part.set_cell(parent.cell)
    for parent_atom in parent:
        found = False
        for part_atom in part:
            if np.allclose(parent_atom.position, part_atom.position, atol=1e-4):
                found = True
                break
        if not found:
            tmp_second_part.append(parent_atom)
    return tmp_second_part

def get_rand_xyz(cell):
    """
    Returns a random xyz coordinates in appropriate ranges depending on the cell size.

    Args:
        cell (ase.Cell): The cell object.

    Returns:
        np.array: The random xyz coordinates.
    """
    array = cell.cellpar()
    a = array[0]  # wymiar w kierunku osi x
    b = array[1]
    c = array[2]  # wymiar w kierunku osi z
    b_y = b * np.sqrt(3) / 2  # skladowa wektora b w kierunku y

    y0 = 0  #zakres dozwolonych wartosci y
    y1 = b_y
    tmp_y = random.uniform(y0, y1)
    x_prim = tmp_y * np.sqrt(3) / 3
    x0 = 0 - x_prim  #zakres dozwolonych wartosci x
    x1 = a - x_prim
    tmp_x = random.uniform(x0, x1)
    tmp_z = c / 2 #wartosc z

    return np.array([tmp_x, tmp_y, tmp_z])

def get_image_positions(cell, atom):
    """
        Returns the positions of atomic images based on the size of the cell.

        Args:
            cell (ase.Cell): The cell object.
            atom (ase.Atom): The atom object.

        Returns:
            np.array: The positions of atom images based on the size of the cell.
    """
    array = cell.cellpar()
    a = array[0]  # Dimension in the x-axis direction
    b = array[1]
    b_y = b * np.sqrt(3) / 2  # Component b in the y direction
    b_x = b / 2  # Component b in the x direction

    tmp_image_positions = [atom.position,
                           atom.position + np.array([a, 0., 0.]),
                           atom.position + np.array([2 * a, 0., 0.]),
                           atom.position + np.array([-b_x, b_y, 0.]),
                           atom.position + np.array([a - b_x, b_y, 0.]),
                           atom.position + np.array([2 * a - b_x, b_y, 0.]),
                           atom.position + np.array([-2 * b_x, 2 * b_y, 0.]),
                           atom.position + np.array([a - 2 * b_x, 2 * b_y, 0.]),
                           atom.position + np.array([2 * a - 2 * b_x, 2 * b_y, 0.])]
    return tmp_image_positions
