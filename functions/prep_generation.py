#prep_generation.py
#funkcja generujaca nowe pokolenie na podstawie populacji podanej w pliku .traj

import os, random
from pathlib import Path
import numpy as np
from ase import io
from ase.io import Trajectory
from math import ceil
from functions.sort_population import sort_population
from functions.gen_rand_struct import gen_rand_struct
from functions.mutation import mutation
from functions.crossover import crossover
from functions.gen_energy_file import gen_energy_file


def prep_generation(pop_filename, pop_size, n_best, n_child, n_mut,
                    struct_filename, size, n_atoms, n_change, atom_symbol,
                    calc, mag_moment, label, new_pop_name):

    #wczytanie populacji poczatkowej z pliku .traj
    previous_pop = sort_population(Trajectory(pop_filename, 'r'))
    better_part = previous_pop[:ceil(pop_size/2)] #lepiej przystosowana czesc populacji

    folder_path = Path(f'{new_pop_name}')
    folder_path.mkdir(parents=True, exist_ok=True)
    original_directory = os.getcwd()
    os.chdir(folder_path)

    new_pop = Trajectory(f'{new_pop_name}.traj', 'w')
    for i in range(n_best):
        new_pop.write(better_part[i]) #dodanie najlepiej przystosowanych osobnikow z poprzedniej populacji
                                      #bezposrednio do nowego pokolenia

    #utworzenie nowych osobnikow w wyniku krzyzowania
    child_counter = 0
    childes = []
    while len(childes) < n_child:
        child_counter += 1
        tmp_folder_path = Path(f'child{child_counter}')
        tmp_folder_path.mkdir(parents=True, exist_ok=True)
        os.chdir(tmp_folder_path)
        parent1, parent2 = random.sample(better_part, 2)
        child = crossover(parent1, parent2, n_change, f'{original_directory}/{struct_filename}', size)
        moments = [0] * (len(child) - n_atoms) + [mag_moment] * n_atoms
        child.set_initial_magnetic_moments(moments)
        try:
            child.calc = calc
            pot_energy = child.get_potential_energy()
            relaxed_child = io.read(f'{label}.XV')
            relaxed_child.pbc = [True, True, False]
            relaxed_child.info['pot_energy'] = np.round(pot_energy, 4)
            childes.append(relaxed_child)
            io.write(f'relaxed_child{child_counter}.xyz', relaxed_child)
            print(f'Successfully relaxed child{child_counter}.')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Successfully relaxed child{child_counter}.\n')
        except Exception as e:
            print(f'Failed to relax child{child_counter}. Error: {e}')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Failed to relax child{child_counter}. Error: {e}\n')


        os.chdir(original_directory / folder_path)

    for child in childes:
        new_pop.write(child)

    #utworzenie nowych osobnikow w wyniku mutacji
    mut_counter = 0
    mut_structures = []
    while len(mut_structures) < n_mut:
        mut_counter += 1
        tmp_folder_path = Path(f'mut{mut_counter}')
        tmp_folder_path.mkdir(parents=True, exist_ok=True)
        os.chdir(tmp_folder_path)
        mut_struct = mutation(random.choice(better_part))
        moments = [0] * (len(mut_struct) - n_atoms) + [mag_moment] * n_atoms
        mut_struct.set_initial_magnetic_moments(moments)
        try:
            mut_struct.calc = calc
            pot_energy = mut_struct.get_potential_energy()
            relaxed_mut_struct = io.read(f'{label}.XV')
            relaxed_mut_struct.pbc = [True, True, False]
            relaxed_mut_struct.info['pot_energy'] = np.round(pot_energy, 4)
            mut_structures.append(relaxed_mut_struct)
            io.write(f'relaxed_mut{mut_counter}.xyz', relaxed_mut_struct)
            print(f'Successfully relaxed mut{mut_counter}.')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Successfully relaxed mut{mut_counter}.\n')
        except Exception as e:
            print(f'Failed to relax mut{mut_counter}. Error: {e}')
            with open(f'../log_{new_pop_name}.txt', 'a') as f:
                f.write(f'Failed to relax mut{mut_counter}. Error: {e}\n')

        os.chdir(original_directory / folder_path)

    for mut_struct in mut_structures:
        new_pop.write(mut_struct)

    #utworzenie pozostalych osobnikow poprzez losowanie nowych struktur
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