#mutation.py
#funkcja operatora mutacji, zmieniajaca pozycje jednego, losowo wybranego atomu znajdujacego sie pomiedzy warstwami
#w przekazanym obiekcie klasy Atoms
#mutation(obiekt_atoms)

from ase import Atoms, Atom
import math
import random
import numpy as np
from functions.small_functions import get_r, calc_distance, get_rand_xyz, get_image_positions


def mutation(atoms):
    structure = atoms.copy() #wyjsciowa struktura bedaca kopia tej podanej jako argument funkcji
    atom_symbol = structure[len(structure)-1].symbol #symbol atomy znajdujacego sie pomiedzy warstwami - potrzebny pozniej

    #wymiary komorki elementarnej
    cell = structure.get_cell()
    array = structure.cell.cellpar()
    a = array[0] #wymiar w kierunku osi x
    b = array[1]
    b_y = b * np.sqrt(3) / 2  #skladowa b w kierunku y
    b_x = b / 2  #skladowa b w kierunku x
    c = array[2] #wymiar w kierunku osi z
    c_half = c / 2 #polowa wysokosci komorki

    atom_indexes = [] #lista z indeksami atomow lezacych pomiedzy warstwami

    #dodanie indeksow atomow lezacych pomiedzy warstwami do listy atom_indexes
    for atom in structure:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            atom_indexes.append(atom.index)

    index_to_del = random.choice(atom_indexes) #losowo wybrany ideks atomu do usuniecia
    del structure[index_to_del] #usuniecie atomu ze struktury
    atom_indexes.remove(index_to_del) #usuniecie indeksu usunietego atomu z tablicy

    #tymczasowa struktura przechowujaca pozycje atomow znajdujacych sie pomiedzy warstwami
    #potrzebna do utworzenia obiektu image_structures
    tmp_structure = Atoms()
    tmp_structure.set_cell(cell)
    for atom in structure:
        if math.isclose(atom.position[2], c_half, abs_tol=1E-1):
            tmp_structure.append(atom)

    tol_r = 0.1 #tolerancja uzywana przy sprawdzaniu odleglosci miedzy dodawanymi atomami

    #tymczasowej struktura 3x3 przechowujacej pozycje atomow znajdujacych sie pomiedzy warstwami i ich obrazy
    image_structure = Atoms()
    image_structure.set_cell(cell)
    image_structure = image_structure * [3, 3, 1]

    for atom in tmp_structure: #dodanie atomow na pozycje obrazow na podstawie obiektu tmp_structures
        image_positions = get_image_positions(cell, atom)

        for position in image_positions:
            image_structure.append(Atom(atom_symbol, position))

    should_continue = True #flaga informujaca czy udalo sie wylosowac pozycje atomu, ktora nie koliduje z obrazami atomow
                           #i atomami z warstw w strukturze docelowej
    while should_continue:
        collision_with_image = False #flaga informujaca o kolizji z obrazami atomow pomiedzy warstwami
        random_positions = get_rand_xyz(cell)
        random_atom = Atom(atom_symbol, random_positions) #atom z nowymi, losowymi wspolrzednymi
        middle_positions = random_positions + np.array([a - b_x, b_y, 0.])
        middle_atom = Atom(atom_symbol, middle_positions) #atom z losowymi wspolrzednymi przesuniety na srodek image_structure

        for atom in image_structure: #sprawdzenie odleglosci pomiedzy middle_atom i obrazami
            d = calc_distance(middle_atom, atom)
            if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_image = True
                    break

        if not collision_with_image:
            collision_with_structure = False #flaga informujaca o kolizji z atomami w strukturze docelowej
            for atom in structure: #sprawdzenie odleglosci pomiedzy random_atom i atomami w strukturze docelowej
                d = calc_distance(random_atom, atom)
                if d < (get_r(random_atom.number) + get_r(atom.number) + tol_r):
                    collision_with_structure = True
                    break

            if not collision_with_structure:
                structure.append(random_atom)
                should_continue = False

    return structure
