# ==================================================
# Imports
# ==================================================
from ase.io import Trajectory
from ase import io
from tqdm import tqdm
from itertools import combinations
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions.prep_struct import prep_struct
from functions.small_functions import *
from functions.load_config import load_config

# ==================================================
# Loading script parameters from the input file
# ==================================================
config = load_config('combi_input.txt')

positions_filename = config['positions_filename']
struct_filename = config['struct_filename']
size = config['size']
n_atoms = config['n_atoms']
draw_range = config['draw_range']


# ==================================================
# The main logic of the program
# ==================================================
def main():
    # Loading the positions where atoms can be added to the 1x1 elementary cell
    pos1x1 = io.read(positions_filename)
    # Positions of atoms in a cell with specified dimensions
    new_pos = pos1x1 * [int(size[0]),int(size[2]),1]
    # Creating a file to store unique structures
    traj_out = Trajectory(f'{size}_{n_atoms}a.traj','w')

    # All combinations of selecting positions for n_atoms from the new_pos object
    combi = list(combinations(new_pos, n_atoms))

    unique_structures = []
    unique_atoms = []
    tol_r = 0.1 # Tolerance used when checking the distance between atoms

    for _ in tqdm(range(draw_range)):
        # Random selection of a single combination
        rand = random.choice(combi)
        # Preparation of a two-layer structure to which atoms will be added
        struct = prep_struct(struct_filename, f'{size}')

        # Creation of an expanded structure to check for collisions between atoms and their images in adjacent cells.
        tmp_struct = Atoms()
        tmp_struct.cell = struct.cell
        for atom in rand:
            tmp_struct.append(atom)
        image_struct = tmp_struct * [3, 3, 1]

        # Checking for collisions between atoms and their images in adjacent cells
        collision = False
        for i in range(len(image_struct)):
            for j in range(i + 1, len(image_struct)):
                distance = image_struct.get_distance(i, j)
                if distance < (2 * get_r(image_struct.get_atomic_numbers()[0]) + tol_r):
                    collision = True
                    break

        if not collision:
            # Adding atoms to the area between layers.
            for atom in rand:
                struct.append(atom)

            # Object used for comparing structures
            matcher = StructureMatcher()
            # Mapping an ase.Atoms object to a pymatgen.Structure object
            structure = AseAtomsAdaptor.get_structure(struct)

            # Comparison of the structure with structures already existing in the database
            is_unique = True
            for unique_structure in unique_structures:
                if matcher.fit(structure, unique_structure):
                    is_unique = False
                    break

            # If there is no collision, and the structure does not duplicate those already existing in the database -
            # add the structure to the appropriate lists and save it to a file containing unique structures
            if is_unique:
                unique_structures.append(structure)
                unique_atoms.append(struct)
                traj_out.write(struct)

    print('Final number of unique structures: ', len(traj_out))