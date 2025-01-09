"""
The module contains a function 'sort_population' that sorts according to the structure energy.
"""

def sort_population(population):
    """
    Sorts the population in order from the lowest to the highest energy.

    Args:
        population (ase.io.Trajectory): A Trajectory object with the population.

    Returns:
        ase.io.Trajectory: The sorted population.
    """
    sorted_structures = sorted(population, key=lambda atoms: atoms.info['pot_energy'])
    return sorted_structures
