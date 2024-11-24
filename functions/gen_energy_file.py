#gen_energy_file.py
#funkcja zapisujaca do pliku .txt energie struktur w podanej populacji, a także sume energii, wartosc srednia i minimalna
#gen_energy_file(plik/obiekt_z_populacja, nazwa_pliku_wyjsciowego)

import numpy as np
from functions.sort_population import sort_population


def gen_energy_file(population, out_filename):
    sorted_population = sort_population(population)
    energies = []
    for i, structure in enumerate(sorted_population):
        pot_energy = structure.info['pot_energy']
        energies.append(pot_energy)
        with open(f'{out_filename}.txt', 'a') as file:
            file.write(f'{i+1}\t{np.round(pot_energy,4)}\n')

    with open(f'{out_filename}.txt', 'a') as file:
        file.write('\n')
        file.write(f'Sum: {np.round(sum(energies),4)}\n')
        file.write(f'Mean: {np.round(np.mean(energies),4)}\n')
        file.write(f'Min: {np.round(np.min(energies),4)}\n')
