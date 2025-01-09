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
    pos1x1 = io.read(positions_filename)
    new_pos = pos1x1 * [int(size[0]),int(size[2]),1]
    traj_out = Trajectory(f'{size}_{n_atoms}a.traj','w')

    combi = list(combinations(new_pos, n_atoms))

    unique_structures = []
    unique_atoms = []
    tol_r = 0.1

    for _ in tqdm(range(draw_range)):
        rand = random.choice(combi)
        struct = prep_struct(struct_filename, f'{size}')
        tmp_struct = Atoms()
        tmp_struct.cell = struct.cell
        for atom in rand:
            tmp_struct.append(atom)

        image_struct = tmp_struct * [3, 3, 1]
        collision = False
        for i in range(len(image_struct)):
            for j in range(i + 1, len(image_struct)):
                distance = image_struct.get_distance(i, j)
                if distance < (2 * get_r(image_struct.get_atomic_numbers()[0]) + tol_r):
                    collision = True
                    break

        if not collision:
            for atom in rand:
                struct.append(atom)

            matcher = StructureMatcher()
            structure = AseAtomsAdaptor.get_structure(struct)

            is_unique = True
            for unique_structure in unique_structures:
                if matcher.fit(structure, unique_structure):
                    is_unique = False
                    break

            if is_unique:
                unique_structures.append(structure)
                unique_atoms.append(struct)
                traj_out.write(struct)

    print('Final number of unique structures: ', len(traj_out))