from ase import io
from ase.io import Trajectory
from pathlib import Path
import os
import numpy as np
from functions.sort_population import sort_population
from functions.gen_rand_struct import gen_rand_struct
from functions.gen_energy_file import gen_energy_file


def gen_random_pop(pop_size, struct_filename, size, n_atoms, atom_symbol, calc, mag_moment, label, new_pop_name):

    folder_path = Path(f'{new_pop_name}')
    folder_path.mkdir(parents=True, exist_ok=True)
    original_directory = os.getcwd()
    os.chdir(folder_path)

    new_pop = Trajectory(f'{new_pop_name}.traj', 'w')

    candidates_counter = 0
    while len(new_pop) < pop_size:
        candidates_counter += 1
        tmp_struct = gen_rand_struct(f'{original_directory}/{struct_filename}', size, atom_symbol, n_atoms)
        moments = [0] * (len(tmp_struct) - n_atoms) + [mag_moment] * n_atoms
        tmp_struct.set_initial_magnetic_moments(moments)
        tmp_folder_path = Path(f'cand{candidates_counter}')
        tmp_folder_path.mkdir(parents=True, exist_ok=True)
        os.chdir(tmp_folder_path)
        try:
            tmp_struct.calc = calc
            pot_energy = tmp_struct.get_potential_energy()
            relaxed_struct = io.read(f'{label}.XV')
            relaxed_struct.pbc = [True, True, False]
            relaxed_struct.info['pot_energy'] = np.round(pot_energy, 4)
            io.write(f'relaxed_cand{candidates_counter}.xyz', relaxed_struct)
            new_pop.write(relaxed_struct)
            print(f'Successfully relaxed cand{candidates_counter}.')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Successfully relaxed cand{candidates_counter}.\n')
        except Exception as e:
            print(f'Failed to relax cand{candidates_counter}. Error: {e}')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Failed to relax cand{candidates_counter}. Error: {e}\n')

        os.chdir(original_directory / folder_path)

    tmp_pop = sort_population(Trajectory(f'{new_pop_name}.traj', 'r'))
    gen_energy_file(tmp_pop, f'../energy_{new_pop_name}') #zapisanie energi nowej populacji do pliku

    os.chdir(original_directory)

    #zapisanie nowej populacji do pliku .traj w kolejnosci od najnizszej do najwyzszej energi
    out_pop = Trajectory(f'sorted_{new_pop_name}.traj', 'w')
    for structure in tmp_pop:
        out_pop.write(structure)

    print(f'The {new_pop_name} is complete!')
    with open(f'{folder_path}/log_{new_pop_name}.txt', 'a') as f:
        f.write(f'The {new_pop_name} is complete!\n')
