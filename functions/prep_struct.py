"""
The module contains a function 'prep_struct' that prepare a two-layer structure.
"""

from ase import Atoms
from ase import io

def prep_struct(struct_filename, size):
    """
    Prepares a two-layer structure of a given size.

    Args:
        struct_filename (str): Name of the file containing dichalcogenide structure.
        size (str): Size of the structure (e.g., 4x4).

    Returns:
        ase.Atoms: The generated two-layer structure.
    """
    # Structure dimensions
    n = int(size[0])
    m = int(size[2])

    structure = Atoms(io.read(struct_filename))
    # Setting a new cell dimension in the direction of the z-axis
    tmp_cell = structure.get_cell()
    tmp_cell[2] *= 1.4
    structure.set_cell(tmp_cell)

    # Temporary structure needed to add atoms
    tmp_structure = structure.copy()

    # Adding atoms of the second layer
    for atom in tmp_structure:
        atom_to_add = atom
        atom_to_add.position[2] = tmp_structure.cell[2, 2] - atom.position[2]
        structure.append(atom_to_add)

    return structure * [n, m, 1]
