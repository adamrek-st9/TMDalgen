#program ewolujący zadaną ilość pokoleń na podstawie początkowej populacji z pliku .traj

from functions.prep_generation import prep_generation
from functions.calculator import get_calc

n_generations = 3     #liczba pokolen
pop_filename = 'pop0.traj' #nazwa pliku z populacja poczatkowa (.traj)
struct_filename = 'MoS2.xyz' #nazwa pliku ze struktura dichalkogenka
size = '3x3'           #rozmiar struktur w populacji (np. 3x3)
n_atoms = 3            #ilosc atomow miedzy warstwami
n_change = 1           #ilosc atomow wymienianych miedzy strukturami podczas krzyzowania
atom_symbol = 'Mo'     #symbol atomow miedzy warstwami
mag_moment = 6.0       #poczatkowy moment magnetyczny nadawany atomom miedzy warstwami
label = 'MoS2'         #etykieta nadawana plikom kalkulatora Siesta (np. MoS2)

#ustawienie kalkulatora wykorzystywanego do obliczen
calc = get_calc(label)

#przygotowanie pierwszego pokolenia na podstawie populacji poczatkowej wczytanej z podanego pliku
prep_generation(pop_filename, struct_filename, size, n_atoms, n_change, atom_symbol,
                calc, mag_moment, label, 'pop1')

#przygotowanie kolejnych pokolen
for i in range(n_generations-1):
    prep_generation(f'sorted_pop{i+1}.traj', struct_filename, size, n_atoms, n_change,
                    atom_symbol, calc, mag_moment, label, f'pop{i+2}')