# ==================================================
# Imports
# ==================================================
from functions.gen_random_pop import gen_random_pop
from functions.prep_generation import prep_generation
from functions.calculator import get_calc
from functions.load_config import load_config

# ==================================================
# Loading algorithm parameters from the input file
# ==================================================
config = load_config('input.txt')

n_generations = config['n_generations']
pop_size = config['pop_size']
n_best = config['n_best']
n_child = config['n_child']
n_mut = config['n_mut']
struct_filename = config['struct_filename']
size = config['size']
n_atoms = config['n_atoms']
n_change = config['n_change']
atom_symbol = config['atom_symbol']
mag_moment = config['mag_moment']
label = config['label']

# ==================================================
# The main logic of the program
# ==================================================
def main():
    # Setting the calculator used for calculations
    calc = get_calc(label)

    # Generating the initial population
    gen_random_pop(pop_size, struct_filename, size, n_atoms, atom_symbol, calc, mag_moment, label, 'pop0')

    # Preparing the next generations
    for i in range(n_generations-1):
        prep_generation(f'sorted_pop{i}.traj', pop_size, n_best, n_child, n_mut,
                        struct_filename, size, n_atoms, n_change, atom_symbol,
                        calc, mag_moment, label, f'pop{i+1}')

if __name__ == '__main__':
    main()