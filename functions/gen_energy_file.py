"""
The module contains the function 'gen_energy_file', which generates a file containing the energy values
of structures in the new population.
"""

import numpy as np
from functions.sort_population import sort_population

def gen_energy_file(population, out_filename):
    """
    Generates a file containing the energy values of structures in the new population,
    sorted from lowest to highest energy. The file also includes the total sum, average,
    and the lowest energy value at the end.

    Args:
        population (ase.io.Trajectory): A Trajectory object with the population.
        out_filename (str): The name of the output file.

    Returns:
        None: The function does not return a value.
    """
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
